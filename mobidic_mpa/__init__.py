#!/usr/bin/env python3
#
# Copyright (C) 2021
#

__author__ = 'Mobidic'
__authors__ = ['Henri Pegeot','Kevin Yauy','Charles Van Goethem','David Baux']
__copyright__ = 'Copyright (C) 2019'
__license__ = 'Academic License Agreement'
__version__ = '1.2.1a'
__email__ = 'c-vangoethem@chu-montpellier.fr'
__status__ = 'dev'

################################################################################
#
# IMPORT
#
################################################################################
import vcf        # read vcf => PyVCF :https://pyvcf.readthedocs.io/en/latest/
import sys        # system command
import os         # os command
import re         # regex
import argparse   # for options
import logging    # logging messages
import subprocess # launch subprocess
import collections
import json
from operator import lt, le, eq, ne, ge, gt

################################################################################
#
# FUNCTIONS
#
################################################################################

############################################################
def getSoftwarePath(software, expected_folder):
    """
    @summary: Returns the path to the software from the expected_folder if it is present or from the PATH environment variable.
    @param software: [str] Name of the software.
    @param expected_folder: [str] The folder where the software it is supposed to be in mSINGS.
    @return: [str] The path of the software.
    """
    path = os.path.join(expected_folder, software)  # Expected path in mSINGS directory
    if not os.path.exists(path):
        path = wich(software)  # Search in PATH
        if path is None:
            raise Exception("The software {} cannot be found in environment.".format(software))
    return path
########################################

########################################
def wich(software):
    """
    @summary: Returns the path to the software from the PATH environment variable.
    @param software: [str] Name of the software.
    @return: [str/None] The path of the software or None if it is not found.
    """
    soft_path = None
    PATH = os.environ.get('PATH')
    for current_folder in reversed(PATH.split(os.pathsep)):  # Reverse PATh elements to kept only the first valid folder
        eval_path = os.path.join(current_folder, software)
        if os.path.exists(eval_path):
            soft_path = eval_path
    return soft_path
############################################################

############################################################
def check_annotation(vcf_infos, annot_dict):
    """
    @summary: Chek if vcf have all annotations required from config file
    @param vcf_infos: [vcf.reader.infos] All fields info from the VCF header
    @param annot_dict: [annot_dict] Dict of fields search from config file
    @return: [None|list] return None or list of fields missing
    """
    error = None

    vcf_keys_needed = [annot_dict[k]["vcf"] for k in annot_dict]

    for k in vcf_keys_needed:
        if k not in vcf_infos:
            if error is None:
                error = list()
            error.append(k)

    return error
########################################

########################################
def check_split_variants(record):
    """
    @summary: Chek if vcf record is correctly splited
    @param record: [vcf.model._record] Current record of the VCF
    @return: [None|list] return None or specific error
    """
    error = None

    if (len(str(record.REF).split(',')) > 1):
        if error is None:
            error = list()
        error.append("Multiple references variant into VCF.")

    if (len(record.ALT) > 1):
        if error is None:
            error = list()
        error.append("Multi allelic variant into VCF.")

    return error
############################################################

############################################################
# Functions to get value into vcf
def cmp(info, values, opt, record=None):
    """
    @summary: Comparator function used in config file
    @param infos: [vcf.reader.infos] All fields info from the VCF header
    @param values: [dict] Dict of values used to return a result
    @param opt: [opt] Dict of options used for this function
    @param record: [vcf.model._record] Current record of the VCF (or None)
    @return: [None] return None or value defined by config file
    """

    reversed = {
        "lt" : False,
        "le" : False,
        "eq" : False,
        "ne" : False,
        "ge" : True,
        "gt" : True
    }

    try:
        od = collections.OrderedDict(sorted(values.items(), reverse=reversed[opt["op"]]))
    except KeyError:
        log.error(
            "Operator 'op' not recognize or missing in function 'cmp'. "
            "Authorized operator : 'lt', 'le', 'eq', 'ne', 'ge', 'gt'. "
            "Please correct your config file or report an issue on github "
            "(https://github.com/mobidic/MPA/issues)"
        )
        sys.exit()
    except AttributeError:
        log.error(
            "Function 'cmp' needs to be associated with 'values'. "
            "Please correct your config file or report an issue on github "
            "(https://github.com/mobidic/MPA/issues)"
        )
        sys.exit()
    except:
        log.error(
            "Unknown error with 'cmp' function. "
            "Please control your config file then report an issue on github "
            "(https://github.com/mobidic/MPA/issues)"
        )
        sys.exit()

    if info is not None:
        for elt in od:
            if eval(opt["op"])(float(info), float(elt)):
                return od[elt]
        return 0
    else:
        return None
########################################

########################################
# Is Indel Surroundig Splice Site
def is_ISSS(info, values, opt, record):
    """
    @summary: Define if indels surrounding a connonic splice site (2pb by default by annovar)
    @param infos: [vcf.reader.infos] All fields info from the VCF header
    @param values: [dict] Dict of values used to return a result
    @param opt: [opt] Dict of options used for this function
    @param record: [vcf.model._record] Current record of the VCF (or None)
    @return: [None|str] return None or "High"
    """

    try:
        result = (values[info] if (values[info] and record.is_indel) else None)
    except KeyError:
        log.error(
            "Field not recognize or missing in vcf in function 'is_ISS'. "
            "Please control your config file and your VCF or report an issue on github "
            "(https://github.com/mobidic/MPA/issues)"
        )
        sys.exit()

    return result
########################################

########################################
def get_first_value(info, values=None, opt=None, record=None):
    """
    @summary: Get the first value from an info field into a vcf or applied a fct to convert this value
    @param infos: [vcf.reader.infos] All fields info from the VCF header
    @param values: [dict] Dict of values used to return a result
    @param opt: [opt] Dict of options used for this function
    @param record: [vcf.model._record] Current record of the VCF (or None)
    @return: [None|str] return None or raw value or converted value from this field
    """

    if opt is None:
        if values is None or info[0] is None:
            return info[0]
        return values[info[0]]

    result = eval(opt["fct"])(info[0], values, opt, record)

    if result is None:
        return opt["default"]
    else:
        return result
########################################

########################################
def get_spliceAI(info, values=None, opt=None, record=None):
    """
    @summary: Get the status of spliceAI score
    @param info: [vcf.reader.infos] All fields info from the VCF header
    @param values: [dict] Dict of values used to return a result
    @param opt: [opt] Dict of options used for this function
    @param record: [vcf.model._record] Current record of the VCF (or None)
    @return: [None|str] return None or status according to the config file
    """
    if info[0] is None:
        return None
    else:
        spliceAI_annot = dict()
        spliceAI_split = info[0].split("\\x3b")
        for annot in spliceAI_split:
            annot_split = annot.split("\\x3d")
            try:
                spliceAI_annot[annot_split[0]] = annot_split[1]
            except IndexError:
                spliceAI_annot[annot_split[0]] = None

        try:
            DS_AG = float(spliceAI_annot["DS_AG"])
            DS_AL = float(spliceAI_annot["DS_AL"])
            DS_DG = float(spliceAI_annot["DS_DG"])
            DS_DL = float(spliceAI_annot["DS_DL"])
        except KeyError:
            log.error(
                "Key used by spliceAI not recognize in annotation. "
                "Key allowed for the resaerch of splicing 'DS_AG', 'DS_AL', 'DS_DG', 'DS_DL'. "
                "Please control your VCF and report an issue on github "
                "(https://github.com/mobidic/MPA/issues)"
            )
            sys.exit()
        except ValueError:
            log.error(
                "The spliceAI score values are not <float>. "
                "Please control your VCF and spliceAI annotation or report an issue on github "
                "(https://github.com/mobidic/MPA/issues)"
            )
            sys.exit()

        maxi = max(DS_AG, DS_AL, DS_DG, DS_DL)

        result = eval(opt["fct"])(maxi, values, opt)
        if result is None:
            return opt["default"]
        else:
            return result
        return result
########################################

########################################
def get_all_infos(record, annot_dict):
    """
    @summary: Get all scores from config file
    @param record: [vcf.model._record] Current record of the VCF (or None)
    @param annot_dict: [annot_dict] Dict of fields search from config file
    @return: [dict] return the dict containing all usefull scores
    """
    r_info = record.INFO
    score_utils = dict()

    for elt in annot_dict:
        # Initialise value to get info in VCF
        vcf_key = annot_dict[elt]["vcf"]
        values = annot_dict[elt]["values"] if "values" in annot_dict[elt] else None
        opt    = annot_dict[elt]["opt"] if "opt" in annot_dict[elt] else None

        # Initialise dict of results
        type_score    = annot_dict[elt]["type"]
        if type_score not in score_utils:
            score_utils[type_score] = list()

        result = eval(annot_dict[elt]["fct"])(r_info[vcf_key], values, opt, record)
        log.debug("{} -> {}".format(elt,r_info[vcf_key]))

        score_utils[type_score].append((elt, result))

    return score_utils
############################################################
# Process scores
def compute_scores(scores, annotation_dict):
    """
    @summary: Compute all score following their weight and penality
    @param record: [scores] Dict of scores associated to their type
    @param annot_dict: [annot_dict] Dict of fields search from config file
    @return: [list] return a list of tuples containing computes scores associated to their type
    """
    compute_scores = list()

    for type in scores:
        # Initialise sum and available counter
        sum_score = 0
        available = 0
        msg = ""
        d_weight = 1
        d_penality = -0.1

        for key, value in scores[type]:
            try:
                weight = annotation_dict[key]["weight"]
            except KeyError:
                weight = d_weight
                log.warning("No weight score for '{}' (default {})".format(key, d_weight))

            try:
                penality = annotation_dict[key]["penality"]
            except KeyError:
                log.warning("No penality score for '{}' (default {})".format(vcf_key, d_penality))

            if value is not None:
                available += annotation_dict[key]["weight"]
                score_weighted = value * annotation_dict[key]["weight"]
            else:
                score_weighted = annotation_dict[key]["penality"]
                weight = 0

            sum_score += score_weighted

            if score_weighted <= 0:
                msg_tmp = "\t- {0:.<18s} {1:<.2f} <= {2:s}".format(key,score_weighted,str(value))
            else:
                msg_tmp = "\t- {0:.<18s} {1:<.2f} <= {2:<.2f} * {3:<.2f}".format(key,score_weighted,value,weight)
            msg = "{}\n{}".format(msg, msg_tmp)

        if available == 0:
            s_adjusted = None
            msg_tmp = "{0:.<10s} {1:s}".format(type, str(s_adjusted))
        else:
            s_adjusted = (sum_score/available) if sum_score/available > 0 else 0
            msg_tmp = "{0:.<10s} {1:<.2f} <= {2:<.2f}/{3:<.2f}".format(type, s_adjusted, sum_score, available)
        msg="{}{}".format(msg_tmp, msg)
        log.debug(msg)

        compute_scores.append((type, s_adjusted))

    return compute_scores


############################################################

################################################################################
#
# PROCESS
#
################################################################################
def main(args, logger):
    """
    @summary: Launch annotation with MPA score on a vcf.
    @param args: [Namespace] The namespace extract from the script arguments.
    @param log: [Logger] The logger of the script.
    """
    global log
    log = logger

    # Get arguments
    input     = args.input
    conf_file = args.configuration
    output    = args.output

    ########################################
    # Create new info fields to add into the vcf output
    # TODO: improve this ! already existing on pyVCF
    _Info = collections.namedtuple(
        'Info', [
            'id',
            'num',
            'type',
            'desc',
            'source',
            'version'
        ]
    )
    info_MPA_adjusted = _Info(
        "MPA_adjusted",
        ".",
        "String",
        "MPA_adjusted : normalize MPA missense score from 0 to 10",
        "MPA",
        __version__
    )
    info_MPA_available = _Info(
        "MPA_available",
        ".",
        "String",
        "MPA_available : number of missense tools annotation available for this variant",
        "MPA",
        __version__
    )
    info_MPA_deleterious = _Info(
        "MPA_deleterious",
        ".",
        "String",
        "MPA_deleterious : number of missense tools that annotate this variant pathogenic",
        "MPA",
        __version__
    )
    info_MPA_final_score = _Info(
        "MPA_final_score",
        ".",
        "String",
        "MPA_final_score : unique score that take into account curated database, biological assumptions, splicing predictions and the sum of various predictors for missense alterations. Annotations are made for exonic and splicing variants up to +300nt.",
        "MPA",
        __version__
    )
    info_MPA_impact = _Info(
        "MPA_impact",
        ".",
        "String",
        "MPA_impact : pathogenic predictions (clinvar_pathogenicity, splice_impact, stop and frameshift_impact)",
        "MPA",
        __version__
    )
    info_MPA_ranking = _Info(
        "MPA_ranking",
        ".",
        "String",
        "MPA_ranking : prioritize variants with ranks from 1 to 10",
        "MPA",
        __version__
    )
    #
    ########################################

    ########################################
    # Get configuration file
    try:
        with open(conf_file) as json_data:
            conf_dict = json.load(json_data)
    except (IsADirectoryError, FileNotFoundError):
        log.error(
            "No such file '{}'. "
            "Please provide a valid file or report an issue on github "
            "(https://github.com/mobidic/MPA/issues)".format(conf_file)
        )
        sys.exit()
    except json.decoder.JSONDecodeError:
        log.error(
            "File '{}' seems to be not a valid JSON file. "
            "Please provide a valid file or report an issue on github "
            "(https://github.com/mobidic/MPA/issues)".format(conf_file)
        )
        sys.exit()
    except:
        log.error(
            "Unexptected error : '{}'. "
            "Please report an issue on github "
            "(https://github.com/mobidic/MPA/issues)".format(sys.exc_info()[0])
        )
        sys.exit()
    #
    ########################################

    ########################################
    # Read the vcf
    try:
        num_lines = sum(1 for line in open(input) if not line.startswith("#") )
        toolbar_width = 50
        keep = num_lines/toolbar_width
        f = open(input, 'r')
    except (IOError):
        log.error(
            "No such file '{}' (IOError). "
            "Please provide a valid file or report an issue on github "
            "(https://github.com/mobidic/MPA/issues)".format(input)
        )
        sys.exit()
    except:
        log.error(
            "Unexptected error : '{}'. "
            "Please report an issue on github "
            "(https://github.com/mobidic/MPA/issues)".format(sys.exc_info()[0])
        )
        sys.exit()
    else:
        with f:
            log.info("Read VCF")
            vcf_reader = vcf.Reader(f)

            # Check the annotation into the fields info in the header of the VCF
            log.info("Check vcf annotations")
            error = check_annotation(vcf_reader.infos, conf_dict["annotations"])
            if error is not None:
                for e in error:
                    log.error(
                        "Expected annotation with info field : '{}'. "
                        "Please refer to the documentation to use Annovar "
                        "or report an issue on github "
                        "(https://github.com/mobidic/MPA/issues)".format(e)
                    )
                sys.exit()

            # TODO: improve this
            vcf_reader.infos.update({'MPA_adjusted':info_MPA_adjusted})
            vcf_reader.infos.update({'MPA_available':info_MPA_available})
            vcf_reader.infos.update({'MPA_deleterious':info_MPA_deleterious})
            vcf_reader.infos.update({'MPA_final_score':info_MPA_final_score})
            vcf_reader.infos.update({'MPA_impact':info_MPA_impact})
            vcf_reader.infos.update({'MPA_ranking':info_MPA_ranking})

            # Initialise output VCF
            vcf_writer = vcf.Writer(open(output, 'w'), vcf_reader)

            toolbar_width = 50
            log.info("Read each variants")
            sys.stdout.write("[%s]" % (" " * toolbar_width))
            sys.stdout.flush()
            sys.stdout.write("\b" * (toolbar_width+1))
            cnt_var = 0
            cnt_write = 0
            for record in vcf_reader:
                cnt_var += 1
                if cnt_var >= keep :
                    cnt_var = 0
                    t = keep
                    while t < 1:
                        sys.stdout.write("-")
                        sys.stdout.flush()
                        cnt_write += 1
                        t += keep
                    if(cnt_write < toolbar_width):
                        cnt_write += 1
                        sys.stdout.write("-")
                        sys.stdout.flush()

                # Confirm that the variant is correctly splitted
                error = check_split_variants(record)
                if error is not None:
                    for e in error:
                        log.error(
                            "Variant '{}' seems not splited : '{}'. "
                            "Please refer to the documentation to split your VCF "
                            "or report an issue on github "
                            "(https://github.com/mobidic/MPA/issues)".format(e)
                        )
                    sys.exit()

                # Get all scores
                scores = get_all_infos(record, conf_dict["annotations"])
                log.debug(scores)
                c_scores = compute_scores(scores, conf_dict["annotations"])
                log.debug(c_scores)
                max_score = ("",0)
                m = 0
                t = 0
                for sc in c_scores:
                    if sc[1] is not None:
                        if (max_score[1] < sc[1]):
                            max_score=sc

                        m += sc[1] * conf_dict["type"][sc[0]]["weight"]
                        t += conf_dict["type"][sc[0]]["weight"]
                    else:
                        m += conf_dict["type"][sc[0]]["penality"]
                        t += 1

                record.INFO['MPA_impact'] = max_score[0]
                record.INFO['MPA_final_score'] = max_score[1]
                record.INFO['MPA_deleterious'] = c_scores
                try:
                    record.INFO['MPA_ranking'] = float(m)/float(t) if float(m)/float(t) >=0 else 0
                except ZeroDivisionError:
                    record.INFO['MPA_ranking'] = None

                vcf_writer.write_record(record)
            vcf_writer.close()
    while(cnt_write < toolbar_width):
        cnt_write += 1
        sys.stdout.write("-")
        sys.stdout.flush()

    sys.stdout.write("]\n")
