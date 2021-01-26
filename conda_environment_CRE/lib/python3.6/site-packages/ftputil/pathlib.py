# Copyright (C) 2002-2019, Stefan Schwarzer <sschwarzer@sschwarzer.net>
# and ftputil contributors (see `doc/contributors.txt`)
# See the file LICENSE for licensing terms.

import os


__all__ = ["FTPPath"]


# class FTPPath(os.PathLike):
#
#     # Needed methods
#     def __init__(self, path):
#         # `path` can be `str`, `bytes` or `PathLike`
#         parts
#         drive
#         root
#         anchor
#         parents
#         parent
#         name
#         suffix
#         suffixes
#         stem
#
#     `PurePath` methods
#     __truediv__
#     as_posix
#     as_uri
#     is_absolute
#     is_reserved
#     joinpath
#     match
#     relative_to
#     with_name
#     with_suffix
#
#     non-`PurePath` methods
#     cwd
#     home
#     stat
#     chmod
#     exists
#     expanduser
#     glob
#     group
#     is_dir
#     is_file
#     is_mount
#     is_symlink
#     is_socket
#     is_fifo
#     is_block_device
#     is_char_device
#     iterdir
#     lchmod
#     lstat
#     mkdir
#     open
#     owner
#     read_bytes
#     read_text
#     rename
#     replace
#     resolve
#     rglob
#     rmdir
#     samefile
#     stat
#     symlink_to
#     touch
#     unlink
#     link_to
#     write_bytes
#     write_text
