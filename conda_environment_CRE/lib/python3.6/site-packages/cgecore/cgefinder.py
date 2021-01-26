#!/usr/bin/env python3
import subprocess
import re
import os.path
import sys

# TODO import blaster and make blaster function in CGEFinder
# from cge.blaster.blaster import Blaster


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


class FinderResult():
    def __init__(self, results, align_sbjct=None, align_query=None,
                 align_homo=None):
        self.results = results  # Results
        self.gene_align_query = align_query  # Sequence alignment lines
        self.gene_align_homo = align_homo  # Sequence alignment homolog string
        self.gene_align_sbjct = align_sbjct  # Sequence alignment allele string


class CGEFinder():

    # Variables used by methods to distinguish results created by different
    # methods.
    TYPE_BLAST = "blast"
    TYPE_KMA = "kma"

    @staticmethod
    def kma(inputfile_1, out_path, databases, db_path_kma, min_cov=0.6,
            threshold=0.9, kma_path="cge/kma/kma", sample_name="",
            inputfile_2=None, kma_mrs=None, kma_gapopen=None,
            kma_gapextend=None, kma_penalty=None, kma_reward=None,
            kma_apm=None, kma_memmode=False, kma_nanopore=False, debug=False,
            kma_add_args=None, kma_cge=False, kma_1t1=False):
        """
           TODO: Result storage - Too complex. Not effective.
           Before changing code: Check downstream dependencies of results
                                 dicts.
           Currently the code stores results in four different dicts. This can
           be reduced to just one.
           The main dict that stores all results currently stores one result pr
           hit, and distiguishes hits to the same gene by adding an increasing
           integer (obtained by getting the length of the internal gene dict).
           The main dict also stores an internal dict for each 'gene'
           containing an empty dict for each hit.
           Solution: Create a main<dict> -> gene<list> -> hit<dict> design, and
           remove all the references main -> hit. This reduces redundancy and
           the need for manipulating gene names with an incremental integer.

           TODO: Method too many responsibilities
           Solution: Create KMA class, create additional functions.

           Original comment:
           "I expect that there will only be one hit pr gene, but if there are
           more, I assume that the sequence of the hits are the same in the res
           file and the aln file."
           Not sure if this holds for the current code.
        """
        threshold = threshold * 100
        min_cov = min_cov * 100

        kma_results = dict()
        kma_results["excluded"] = dict()

        if(sample_name):
           sample_name = "_" + sample_name

        # Initiate output dicts.
        gene_align_sbjct = {}
        gene_align_query = {}
        gene_align_homo = {}

        for db in databases:
            kma_db = db_path_kma + "/" + db
            kma_outfile = out_path + "/kma_" + db + sample_name
            kma_cmd = ("%s -t_db %s -o %s -e 1.0" % (kma_path,
                       kma_db, kma_outfile))
            if(inputfile_2 is not None):
                kma_cmd += " -ipe " + inputfile_1 + " " + inputfile_2
            else:
                kma_cmd += " -i " + inputfile_1
            if(kma_mrs is not None):
                kma_cmd += " -mrs " + str(kma_mrs)
            if(kma_gapopen is not None):
                kma_cmd += " -gapopen " + str(kma_gapopen)
            if(kma_gapextend is not None):
                kma_cmd += " -gapextend " + str(kma_gapextend)
            if(kma_penalty is not None):
                kma_cmd += " -penalty " + str(kma_penalty)
            if(kma_reward is not None):
                kma_cmd += " -reward " + str(kma_reward)
            if(kma_apm is not None):
                kma_cmd += " -apm " + kma_apm
            if(kma_cge):
                kma_cmd += " -cge "
            if(kma_1t1):
                kma_cmd += " -1t1 "
            if (kma_memmode):
                kma_cmd += " -mem_mode "
            if (kma_nanopore):
                kma_cmd += " -bcNano "
                kma_cmd += " -mp 20 "
            if (kma_add_args is not None):
                kma_cmd += " " + kma_add_args + " "

            # kma output files
            align_filename = kma_outfile + ".aln"
            res_filename = kma_outfile + ".res"

            # If .res file exists then skip mapping
            if(os.path.isfile(res_filename)
               and os.access(res_filename, os.R_OK)):
                eprint("Found " + res_filename + " skipping DB.")
            else:
                # Call KMA
                if(debug):
                    eprint("KMA cmd: " + kma_cmd)

                process = subprocess.Popen(kma_cmd, shell=True,
                                           stdout=subprocess.PIPE,
                                           stderr=subprocess.PIPE)
                out, err = process.communicate()

            kma_results[db] = 'No hit found'

            # Open res file
            try:
                res_file = open(res_filename, "r")
                header = res_file.readline()
            except IOError as error:
                sys.exit("Error: KMA did not run as expected.\n"
                         "KMA finished with the following response:"
                         "\n{}\n{}".format(out.decode("utf-8"),
                                           err.decode("utf-8")))

            gene_res_count = {}

            for line in res_file:

                if kma_results[db] == 'No hit found':
                    kma_results[db] = dict()

                data = [data.strip() for data in line.split("\t")]
                gene = data[0]

                gene_count = gene_res_count.get(gene, 0)
                gene_count = gene_count + 1
                gene_res_count[gene] = gene_count
                if(gene_count == 1):
                    hit = gene
                else:
                    hit = "{}_{}".format(gene, gene_count)

                sbjct_len = int(data[3])
                sbjct_ident = float(data[4])
                coverage = float(data[5])
                depth = float(data[-3])
                q_value = float(data[-2])
                p_value = float(data[-1])

                exclude_reasons = []

                if(coverage < min_cov or sbjct_ident < threshold):
                    exclude_reasons.append(coverage)
                    exclude_reasons.append(sbjct_ident)

                if(exclude_reasons):
                    kma_results["excluded"][hit] = exclude_reasons

                kma_results[db][hit] = dict()
                kma_results[db][hit]['sbjct_length'] = sbjct_len
                kma_results[db][hit]["perc_coverage"] = coverage
                kma_results[db][hit]["sbjct_string"] = []
                kma_results[db][hit]["query_string"] = []
                kma_results[db][hit]["homo_string"] = []
                kma_results[db][hit]["sbjct_header"] = gene
                kma_results[db][hit]["perc_ident"] = sbjct_ident
                kma_results[db][hit]["query_start"] = "NA"
                kma_results[db][hit]["query_end"] = "NA"
                kma_results[db][hit]["contig_name"] = "NA"
                kma_results[db][hit]["HSP_length"] = ""
                kma_results[db][hit]["cal_score"] = q_value
                kma_results[db][hit]["depth"] = depth
                kma_results[db][hit]["p_value"] = p_value
            res_file.close()

            if kma_results[db] == 'No hit found':
                continue

            # Open align file
            with open(align_filename, "r") as align_file:
                gene_aln_count = {}
                gene = ""
                # Parse through alignments
                for line in align_file:
                    # Skip empty lines
                    if(not line.strip()):
                        continue

                    # Check when a new gene alignment start
                    if line.startswith("#"):
                        gene = line[1:].strip()

                        gene_count = gene_aln_count.get(gene, 0)
                        gene_aln_count[gene] = gene_count + 1
                    else:
                        if(gene_aln_count[gene] == 1):
                            hit = gene
                        else:
                            hit = "{}_{}".format(gene, gene_aln_count[gene])

                        if hit in kma_results[db]:
                            line_data = line.split("\t")[-1].strip()
                            if line.startswith("template"):
                                kma_results[db][hit]["sbjct_string"] += (
                                    [line_data])
                            elif line.startswith("query"):
                                kma_results[db][hit]["query_string"] += (
                                    [line_data])
                            else:
                                kma_results[db][hit]["homo_string"] += (
                                    [line_data])
                        else:
                            print(hit + " not in results: ", kma_results)

            # concatinate all sequence lists and find subject start
            # and subject end

            gene_align_sbjct[db] = {}
            gene_align_query[db] = {}
            gene_align_homo[db] = {}

            for hit in kma_results[db]:
                align_sbjct = "".join(kma_results[db][hit]['sbjct_string'])
                align_query = "".join(kma_results[db][hit]['query_string'])
                align_homo = "".join(kma_results[db][hit]['homo_string'])

                # Extract only aligned sequences
                start = re.search("^-*(\w+)", align_query).start(1)
                end = re.search("\w+(-*)$", align_query).start(1)

                kma_results[db][hit]['sbjct_string'] = align_sbjct[start:end]
                kma_results[db][hit]['query_string'] = align_query[start:end]
                kma_results[db][hit]['homo_string'] = align_homo[start:end]

                # Save align start and stop positions relative to
                # subject sequence
                kma_results[db][hit]['sbjct_start'] = start + 1
                kma_results[db][hit]["sbjct_end"] = end + 1
                kma_results[db][hit]["HSP_length"] = end - start

                # Count gaps in the alignment
                kma_results[db][hit]["gaps"] = (
                    kma_results[db][hit]['sbjct_string'].count("-")
                    + kma_results[db][hit]['query_string'].count("-"))

                # Save sequences covering the entire subject sequence
                # in seperate variables
                gene_align_sbjct[db][hit] = align_sbjct
                gene_align_query[db][hit] = align_query
                gene_align_homo[db][hit] = align_homo

        return FinderResult(kma_results, gene_align_sbjct, gene_align_query,
                            gene_align_homo)
