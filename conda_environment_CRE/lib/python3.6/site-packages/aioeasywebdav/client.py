import os
import ssl
import time
import shutil
import logging
import hashlib
import asyncio
import aiohttp
import platform
import threading
from abc import ABCMeta, abstractmethod
from numbers import Number
import xml.etree.cElementTree as xml
from collections import namedtuple
from urllib.parse import urlparse, quote, unquote

from http.client import responses as HTTP_CODES

DOWNLOAD_CHUNK_SIZE_BYTES = 1 * 1024 * 1024
MAX_OPEN_FILES = 512
TEMP_NAME='.part'


class WebdavException(Exception):
    pass


class ConnectionFailed(WebdavException):
    pass


def codestr(code):
    return HTTP_CODES.get(code, 'UNKNOWN')


File = namedtuple('File', ['name', 'size', 'mtime', 'ctime', 'contenttype', 'checksum'])


def prop(elem, name, default=None, ns='{DAV:}'):
    child = elem.find('.//' + ns + name)
    return default if child is None else child.text


def elem2file(elem, baseurl='', basepath=''):
    length = prop(elem, 'getcontentlength', 0)
    return File(
        unquote(prop(elem, 'href').replace(baseurl, '').replace(basepath, '')),
        int(length) if length else 0,
        prop(elem, 'getlastmodified', ''),
        prop(elem, 'creationdate', ''),
        prop(elem, 'getcontenttype', ''),
        # Owncloud Webdav checksum support. Unsure if regular owncloud client creates the checksum, if not refer to:
        # https://github.com/owncloud/core/pull/21997#issuecomment-176972659
        # Silently ignored on other webdav servers.
        prop(elem, 'checksum', None, ns='{http://owncloud.org/ns}'),
    )


class OperationFailed(WebdavException):
    _OPERATIONS = dict(
        HEAD = "get header",
        GET = "download",
        PUT = "upload",
        DELETE = "delete",
        MKCOL = "create directory",
        PROPFIND = "list directory",
        )

    def __init__(self, method, path, expected_code, actual_code):
        self.method = method
        self.path = path
        self.expected_code = expected_code
        self.actual_code = actual_code
        operation_name = self._OPERATIONS[method]
        self.reason = 'Failed to {operation_name} "{path}"'.format(**locals())
        expected_codes = (expected_code,) if isinstance(expected_code, Number) else expected_code
        expected_codes_str = ", ".join('{0} {1}'.format(code, codestr(code)) for code in expected_codes)
        actual_code_str = codestr(actual_code)
        msg = '''\
{self.reason}.
  Operation     :  {method} {path}
  Expected code :  {expected_codes_str}
  Actual code   :  {actual_code} {actual_code_str}'''.format(**locals())
        super(OperationFailed, self).__init__(msg)


class ProgressHandler(object):
    __metaclass__ = ABCMeta

    @abstractmethod
    async def progress_callback(self, total_size, current_size):
        """
        :param int total_size: Eventual total size of the file
        :param int current_size: Current downloaded size of the file
        :return: None
        """
        pass

    @abstractmethod
    async def done_callback(self, error):
        """
        Called when the download task is finished
        :param (str | None) error: Error string if there was an issue or None
        :return: None
        """
        pass

    @abstractmethod
    async def verifying_callback(self):
        """
        Called when the download transfer is complete and about to be verified
        :return: None
        """
        pass


class Client(object):

    def __init__(self, url=None, host=None, port=0, auth=None, username=None, password=None,
                 protocol='http', verify_ssl=True, path=None, cert=None, max_connections=65):
        self.log = logging.getLogger("%s.%s" % (self.__class__.__module__, self.__class__.__name__))
        self._max_connections = max_connections

        self.basepath = ''
        if url:
            self.baseurl = url
            self.basepath = urlparse(url).path
        else:
            if not port:
                port = 443 if protocol == 'https' else 80
            self.baseurl = '{0}://{1}:{2}'.format(protocol, host, port)
            if path:
                self.baseurl = '{0}/{1}'.format(self.baseurl, path)
                self.basepath = path
        self.cwd = '/'

        sslcontext = None
        if cert:
            sslcontext = ssl.create_default_context(cafile=cert)
        conn = aiohttp.TCPConnector(ssl_context=sslcontext, verify_ssl=verify_ssl, limit=self._max_connections)

        if not auth and username and password:
            auth = aiohttp.BasicAuth(username, password=password)

        self.session = aiohttp.ClientSession(connector=conn, auth=auth)

        self._closed = False
        self._rate_ave_period = 2  # approx period in seconds over which the download rate it averaged
        self._download_rates = {}
        self._rate_tracking = {}
        self._rate_calc_future = asyncio.ensure_future(self._rate_calc())

        self.limit_files = asyncio.Semaphore(MAX_OPEN_FILES)

    async def close(self):
        self.session.close()
        self._closed = True
        await self._rate_calc_future

    async def _rate_calc(self):
        prev = time.time() - self._rate_ave_period
        while not self._closed:
            now = time.time()
            for name, downloaded in self._rate_tracking.items():
                # assert isinstance(entry, self._rate_tracking_entry)
                period = now - prev
                # moving = self._download_rates.get(name, 0)
                # moving = moving + downloaded - moving/period
                # self._download_rates[name] = moving
                self._download_rates[name] = downloaded / period
                self._rate_tracking[name] = 0
            prev = now
            await asyncio.sleep(self._rate_ave_period)

    def _rate_notify(self, name, downloaded):
        self._rate_tracking[name] = self._rate_tracking.get(name, 0) + downloaded
        self._rate_tracking[None] = self._rate_tracking.get(None, 0) + downloaded

    def download_rate(self, name = None):
        """
        :param (str or None) name: filename or None for global rate
        :return: download rate in bytes per second
        """
        return self._download_rates.get(name, 0)# / self._rate_ave_period

    async def _send(self, method, path, expected_code, **kwargs):
        url = self._get_url(path)
        response = await self.session.request(method, url, allow_redirects=False, **kwargs)
        if isinstance(expected_code, Number) and response.status != expected_code \
            or not isinstance(expected_code, Number) and response.status not in expected_code:
            status = response.status
            response.close()
            raise OperationFailed(method, path, expected_code, status)
        return response

    def _get_url(self, path):
        path = quote(str(path).strip())
        if path.startswith('/'):
            return self.baseurl + path
        return "".join((self.baseurl, self.cwd, path))

    def cd(self, path):
        path = path.strip()
        if not path:
            return
        stripped_path = '/'.join(part for part in path.split('/') if part) + '/'
        if stripped_path == '/':
            self.cwd = stripped_path
        elif path.startswith('/'):
            self.cwd = '/' + stripped_path
        else:
            self.cwd += stripped_path

    async def mkdir(self, path, safe=False):
        expected_codes = 201 if not safe else (201, 301, 405)
        async with (await self._send('MKCOL', path, expected_codes)):
            pass

    async def mkdirs(self, path):
        dirs = [d for d in path.split('/') if d]
        if not dirs:
            return
        if path.startswith('/'):
            dirs[0] = '/' + dirs[0]
        old_cwd = self.cwd
        try:
            for dir in dirs:
                try:
                    await self.mkdir(dir, safe=True)
                except Exception as e:
                    if e.actual_code == 409:
                        raise
                finally:
                    self.cd(dir)
        finally:
            self.cd(old_cwd)

    async def rmdir(self, path, safe=False):
        path = str(path).rstrip('/') + '/'
        expected_codes = 204 if not safe else (204, 404)
        async with (await self._send('DELETE', path, expected_codes)):
            pass

    async def delete(self, path):
        """
        :param (File | str) path: path from server root or File object to delete
        :return: None
        """
        path = path.name if isinstance(path, File) else path
        async with (await self._send('DELETE', path, 204)):
            pass

    async def upload(self, local_path_or_fileobj, remote_path):
        if isinstance(local_path_or_fileobj, str):
            async with self.limit_files:
                with open(local_path_or_fileobj, 'rb') as f:
                    await self._upload(f, remote_path)
        else:
            await self._upload(local_path_or_fileobj, remote_path)

    async def _upload(self, fileobj, remote_path, **kwargs):
        async with self.limit_files:
            loop = asyncio.get_event_loop()
            local_hash = await loop.run_in_executor(None, self._md5, fileobj)
        headers = {"OC-Checksum" : "MD5:%s" % local_hash}
        async with (await self._send('PUT', remote_path, (200, 201, 204), data=fileobj, headers=headers, **kwargs)):
            pass

    async def _check_existing_download(self, fileobj, remote_file, start):
        overlap = 16
        length = fileobj.seek(0, os.SEEK_END)
        valid = False if length else True
        if length:
            pos = fileobj.seek(-overlap, os.SEEK_CUR) + start
            tail = fileobj.read()
            async with (await self._send(
                    'GET', remote_file.name, (200,206), headers={"Range": "bytes=%d-%d"%(pos,pos+overlap-1)})) as rsp:
                read = await rsp.content.read()
                if read == tail:
                    valid = True
                else:
                    fileobj.seek(0)
        return valid

    async def _download_stream(self, local_path, part_path, remote_file, start, end, progress_callback = None, enabled_event = None):
        finished = False
        pos = existing = 0
        retry = 5
        while not finished and retry:
            retry -= 1
            try:
                if enabled_event:
                    await enabled_event.wait()

                if isinstance(local_path, str):
                    exists = os.path.exists(part_path)
                    mode = 'r+b' if exists else 'w+b'

                    async with self.limit_files:
                        with open(part_path, mode) as fileobj:

                            if exists:
                                await self._check_existing_download(fileobj, remote_file, start)
                            pos = existing = fileobj.tell()
                            if progress_callback:
                                await progress_callback(part_path, pos)
                            req_start = start + pos
                            if req_start >= end:
                                finished = True
                else:
                    req_start = start

                if not isinstance(end, int) or req_start < end:

                    header = {"Range": "bytes=%s-%s" % (req_start, end)}
                    async with (await self._send(
                            'GET', remote_file.name, (200,206), chunked=True, headers=header)) as response:

                        while not enabled_event or enabled_event.is_set():

                            chunk = await response.content.read(DOWNLOAD_CHUNK_SIZE_BYTES)
                            if not chunk:
                                finished = True
                                break

                            self._rate_notify(str(local_path), len(chunk))

                            if isinstance(local_path, str):
                                async with self.limit_files:
                                    with open(part_path, 'r+b') as fileobj:
                                        fileobj.seek(pos)
                                        length = len(chunk)
                                        fileobj.write(chunk)
                                        pos += length
                                        if progress_callback:
                                            await progress_callback(part_path, pos)
                            else:
                                length = len(chunk)
                                local_path.write(chunk)
                                pos += length
                                if progress_callback:
                                    await progress_callback(part_path, pos)


            except Exception as ex:
                await asyncio.sleep(3) # Rate limiting
                if not retry:
                    raise  # somewhere to breakpoint
        return part_path, pos - existing

    @staticmethod
    def join_parts(local_path, partfiles):
        # Join all parts together
        with open(local_path, 'wb') as local_file:
            for part in partfiles:
                with open(part, 'rb') as partfile:
                    shutil.copyfileobj(partfile, local_file)
        for part in partfiles:
            os.unlink(part)

    @staticmethod
    def md5(fname):
        with open(fname, "rb") as f:
            return Client._md5(f)

    @staticmethod
    def _md5(fileobj):
        start = fileobj.tell()
        hash_md5 = hashlib.md5()
        for chunk in iter(lambda: fileobj.read(65535), b""):
            hash_md5.update(chunk)
        fileobj.seek(start)
        return hash_md5.hexdigest()

    async def download(self, remote_file, local_path, progress_handler = None, enabled_event = None):
        """
        :param str local_path: path to write local file
        :param File remote_file: remote file description as returned by ls
        :param (ProgressHandler or None) progress_handler: optional class of callbacks which gets notified of progress
        :param (asyncio.Event or None) enabled_event: optional Event used to pause/resume download
        :return:
        """
        error = 'Unknown'
        responses = []
        cur_lengths = {}
        loop = asyncio.get_event_loop()

        try:
            if isinstance(remote_file, str):
                # Lookup remote file details
                remote_file = (await self.ls(remote_file))[0]

            expected_length = remote_file.size

            chunks = []
            range_details = None

            if isinstance(local_path, str):
                # Download to local filename in multiple part files

                dirname = os.path.dirname(local_path)
                if not dirname:
                    dirname = '.'
                if not os.path.exists(dirname):
                    os.makedirs(dirname)

                chunksize = (10*1000*1000)
                chunks = [(idx, rng, rng + chunksize-1) for idx, rng in enumerate(range(0, remote_file.size, chunksize))]

                if chunks:
                    partname = "%s.{}%s" % (local_path, TEMP_NAME)
                    chunks[-1] = (chunks[-1][0],chunks[-1][1], expected_length)
                    range_details = [(partname.format(str(i)),start,end) for i,start,end in chunks]

            else:  # Assume local_path is a file-like object to stream into
                range_details = [("stream", 0, remote_file.size)]

            if range_details:

                if isinstance(progress_handler, ProgressHandler):
                    async def progress_callback(part_path, cur_length):
                        cur_lengths[part_path] = cur_length
                        await progress_handler.progress_callback(expected_length, sum(cur_lengths.values()))
                else:
                    progress_callback = None

                downloads = [self._download_stream(local_path, part_path, remote_file, start, end, progress_callback, enabled_event)
                             for part_path, start, end in range_details]
                responses = await asyncio.gather(*downloads)

            partfiles = [part for part, dl_size in responses]

            if expected_length == 0:
                error = "Null file reported"
            else:
                if isinstance(local_path, str):  # Check and re-assemble part files
                    local_temp = local_path + TEMP_NAME
                    if sum([os.path.getsize(f) for f in partfiles]) == expected_length and len(partfiles) == len(chunks):
                        error = None
                        if len(chunks) > 1:
                            async with self.limit_files:
                                await loop.run_in_executor(None, self.join_parts, local_temp, partfiles)
                        else:
                            os.rename(partfiles[0], local_temp)

                    if remote_file.checksum:
                        if isinstance(progress_handler, ProgressHandler):
                            await progress_handler.verifying_callback()
                        kind, hash = remote_file.checksum.split(":")
                        local_hash = None
                        if kind == 'MD5':
                            async with self.limit_files:
                                local_hash = await loop.run_in_executor(None, self.md5, local_temp)
                        if local_hash != hash:
                            error = "Invalid Checksum, expected: %s" % str(remote_file.checksum)
                            os.rename(local_temp, local_path+".invalid")

                    if not error:
                        os.rename(local_temp, local_path)

        except Exception as ex:
            tb = __import__('traceback').format_exc()
            error = str(ex)
            if not error:
                error = type(ex)

        finally:
            if isinstance(progress_handler, ProgressHandler):
                await progress_handler.done_callback(error)
            if local_path in self._rate_tracking:
                del self._rate_tracking[local_path]
            if local_path in self._download_rates:
                del self._download_rates[local_path]
        return error

    async def ls(self, remote_path=''):
        """
        :param str remote_path: path relative to the server to list
        :return: [File]
        """
        headers = {'Depth': '1'}
        checksum_prop = """<?xml version="1.0"?>
        <a:propfind xmlns:a="DAV:">
            <a:prop xmlns:oc="http://owncloud.org/ns">
                <oc:checksums/>
                <a:getcontentlength/>
                <a:getlastmodified/>
                <a:creationdate/>
                <a:getcontenttype/>
            </a:prop>
        </a:propfind>"""
        async with (await self._send(
                'PROPFIND', remote_path, (207, 301), headers=headers, data=checksum_prop)) as response:
            # Redirect
            if response.status == 301:
                url = urlparse(response.headers['location'])
                return self.ls(url.path)
            tree = xml.fromstring(await response.read())

        entries = [elem2file(elem, self.baseurl, self.basepath) for elem in tree.findall('{DAV:}response')]
        return entries

    async def exists(self, remote_path):
        async with (await self._send('HEAD', remote_path, (200, 301, 404))) as response:
            ret =  True if response.status != 404 else False
        return ret
