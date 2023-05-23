import sys
import os
import re
import json
import http
import argparse as ap

from Bio import Entrez, SeqIO 
import pycountry_convert as pc
from datetime import datetime as dt

# This script uses NCBI accession number prefixes to remove 'patent' sequence
# such as lab strains and vaccine strains. A list of prefixes can be found
# at https://www.ncbi.nlm.nih.gov/genbank/acc_prefix/

# Patent sequence prefixes from DDBJ
DDBJ_Patents = ["E","BD","DD","DI","DJ","DL","DM","FU","FV", \
    "FW","FZ","GB","HV","HW","HZ","LF","LG", "LV","LX","LY", \
    "LZ","MA","MB","MC","MD","ME","OF","OG","OI","OJ","PA","PB", \
    "PC","PD","PE"]

# Patent sequence prefixes from EMBL
EMBL_Patents = ["A","AX","CQ","CS","FB","GM","GN","HA","HB","HC", \
    "HD","HH","HI","JA","JB","JC","JD", "JE","LP","LQ","MP","MQ","MR","MS"]

# Patent sequene prefixes from Genbank
GB_Patents = ["I","AR","DZ","EA","GC","GP","GV","GX","GY","GZ", \
    "HJ","HK","HL","KH","MI","MM","MO","MV","MX","MY","OO"]

# The total list of prefixes will be the DDBJ, EMBL, and Genbank
# lists combined.
patent_sequence_prefixes = DDBJ_Patents + EMBL_Patents + GB_Patents

# OLD LIST USED - Figured out sequences were missing from this list
"""["A","AX","CQ","CS", "E","FB", \
    "GM","GN","HA","HB","HC","HD","HH","HI","JA","JB","JC","JD", \
    "JE","LP","LQ","MP","MQ","MR","MS","I","AR","DZ","EA","GC", \
    "GP","GV","GX","GY","GZ","HJ","HK","HL","KH","MI","MM","MO", \
    "MV","MX","MY","OO"]"""

# A regex pattern used to grab the prefix from an accession number
accession_regex_pattern = '([A-Z]{1}|[A-Z]{2})_?(\d{6}|\d{5})\.\d'

# In some cases, a sequence will not contain a segment feature
# in the genbank record. If not, we can check the sequence header
# for the segment number (as it may be added here). THis 
description_segment_regex_pattern = '(segment:? (\d))'

# If no RefSeq sequences are returned for a given organism, the following title
# qualifiers will be used to query Nucleotide to determine the average genome length.
# Qualifiers will be searched in order (so first item in list will be use first) and
# will stop if any sequences are returned for that result.
title_search_quantifiers = ["complete genome", "complete"]

# Biopython's Entrez Package sets a cap of 20 records to be obtained by default with a maximum
# of 10,000 possible. This variable stores 10,000 to set all searches to return the maximum and avoid
# missing samples due to this return cap.
retmax = 100000000


def searchNucleotide(query):
    """Uses the Biopython Entrez esearch function to search the NCBI
    nucleotide database using a user supplied query. The esearch function
    returns a handle which the function then parses and returns the parsed entry.

    Parameters
    ----------
    query: str
        The query to be supplied to NCBI

    """

    # Runs the esearch command against the nucleotide database using
    # the provided query
    recordsHandle = Entrez.esearch(db="nucleotide", term=query, retmax=retmax)

    # Parses the handle produced by the search and
    # closes the handle.
    records = Entrez.read(recordsHandle)
    recordsHandle.close()

    # Return the parsed records from the search.
    return records


def getSourceFeature(rec):
    """GenBank records contain multiple 'features' which are annotations 
    describing the sequence. There are typically numerous features in a
    GenBank record, but one important feature is the 'source' feature. This
    contains summary information about the sequence (length, taxID, collection date, etc.).
    Thus, this function searches through the list of features from a sequence and returns the
    source feature.

    Parameters
    ----------
    rec: object (Biopython.SeqRecord)
        A GenBank Biopython SeqRecord object to be parsed
    """

    # Sets the source equal to None (in case no source is found.)
    source = None

    # Loops over all of the features in the record.
    for f in rec.features:
                
        # If the source feature is found.
        if f.type == 'source':
            # Save the source feature in the source variable
            # and break from the loop.
            source = f
            break

    # Return the source feature in dictionary format.
    return source

def getFirstGeneProduct(rec):
    """FUNCTION STILL UNDER DEVELOPMENT
    """

    product = None
    for f in rec.features:
        if 'product' in f.qualifiers.keys():
            product = f.qualifiers['product'][0]

    return product

def getContinent(country):
    """Uses the pycountry-convert library to determine the continent from a country name.

    Parameters
    ----------
    country: str
        A country to convert into a continent.
    """

    # Creates an empty variable to store the continent name.
    continent_name = ""
    
    # The next set of code is wrapped in a try-catch block because
    # pycountry-convert throws an 'KeyError' exception if the country
    # provided does not exist. This has been found to occur with some countries
    # containing accent characters. The NCBI GenBank record will not always include these.
    # However, pycountry-convert may expect them.
    try:
        # Converts the country name to an alpha2 value
        country_alpha2 = pc.country_name_to_country_alpha2(country)
        
        # Converts the country's alpha2 value into its corresponding continent code.
        country_continent_code = pc.country_alpha2_to_continent_code(country_alpha2)

        # Converts the continent code to the continent name
        continent_name = pc.convert_continent_code_to_continent_name(country_continent_code)
    except KeyError as e:
        # If there is an error while converting the country name to continent name,
        # simply print the error and move on. The user can manually enter this metadata
        # after the script completes. THe continent name will be returned as an empty string.
        print(e)

    # Return the continent name.
    return continent_name


def getReportData(s,r):
    """Grabs fields of interest from the source feature of a GenBank record.
    These fields include:
        1. organism taxID
        2. isolate
        3. country of origin
        4. continent
        5. region
        6. isolation source
        7. collection date
    
    Each value in the source feature is called a 'quantifier' and the
    quantifiers provided in the source feature. Thus, the function checks whether
    the values are provided, and returns them if so.

    As well, the function grabs additional fields of interest from a
    reference feature:
    
        1. journal
        2. authors
        3. title

        
    Parameters
    ----------
    s: dict
        The source feature dictionary from a GenBank record

    r: dict
        A reference feature dictionary from a GenBank record.  
    """

    # Creates empty variables to store all of the
    # data fields we want from the source feature.
    taxID = ''
    isolate = ''
    country = ''
    continent = ''
    region = ""
    iso_source = ''
    col_date = ''

    # Check whether the source feature contains the taxonomy database
    # reference quantifier.
    if 'db_xref' in s.qualifiers:
        # If yes, set the value of the taxID variable.
        taxID = s.qualifiers['db_xref'][0].replace('taxon:', "")

    # Check whether the source feature contains the isolate
    # quantifier.
    if 'isolate' in s.qualifiers:
        # If yes, set the value of the isolate variable.
        isolate = s.qualifiers['isolate'][0]
                    
    # Check whether the source feature contains the 
    # country quantifier.
    if 'country' in s.qualifiers:
        # If yes, then we need to first parse the quantifier value. NCBI GenBank files will sometimes
        # contain the country and the region in the same quantifier separated by a semi-colon 
        # (or semicolon + space ': '). Thus we can split the string at this character to separate the
        # country from region information (if contained in the annotation).
        countryString = s.qualifiers['country'][0].replace(": ", ":").split(":")

        # Grab the country from the split country string. Even if there was no ":" character,
        # the split() function will still return a list, and country will still be first.
        country = countryString[0]

        # If the split country quantifier has more than 1 component (after splitting),
        # then the second value must be the region.
        if len(countryString) > 1:
            # Set the value of the region variable, if yes.
            region = countryString[1]
        
        # Convert the country into the continent.
        continent = getContinent(country)
                    
    # Check whether the source feature contains the isolation source
    # quantifier.
    if 'isolation_source' in s.qualifiers:
        # If yes, set the value of the isolation source variable.
        iso_source = s.qualifiers['isolation_source'][0]

    # Check whether the source feature contains the collection source
    # quantifier.
    if 'collection_date' in s.qualifiers:
        # If yes, set the value of the collection date variable.
        col_date = s.qualifiers['collection_date'][0]

    # Grab the journal, authors, and title values from the reference
    # dictionary.
    journal = r.journal
    authors = r.authors
    title = r.title

    # Return the values obtained. If any were not found, an empty string will be returned for that value.
    return taxID, isolate, continent, country, region, iso_source, col_date, authors, title, journal


def getRecSegment(s, recordDesc):
    """Determines whether a record is a segment (of a segmented virus) by
    checking both the source features and the record description. The segment quantifier
    is not a requirement of the source feature. Thus, there could be cases were the record
    is a segment, but it is not marked as so. One way of catching this case is by checking the
    description of the record, as it may include the segment value in the description. If a segment 
    value is found, it is returned. If no segment quantifier or text is found, 'None' is returned. 
    
    Parameters
    ----------
    s: dict
        The source feature dictionary from a GenBank record

    recordDesc: string
        The record description/title
    """
    
    # Sets the segment to None by default.
    segment = None

    # Checks whether the source feature contains
    # a segment quantifier.
    if 'segment' in s.qualifiers:
        # If yes, set the segment to this value
        segment = s.qualifiers['segment'][0]

    # If the segment was not found the source feature qualifiers,
    # see if the description contains the segment using regex.
    if segment == None:
        
        # Performs a regex search on the description using a pattern
        # set to match the format "segment \d".
        r = re.search(description_segment_regex_pattern, recordDesc)
        
        # If there was match, grab the captured segment value and set
        # the segment variable.
        if r:
            segment = r.group(2)

    # Return the segment value.
    return segment


def checkDateFormat(date):
    """Checks the format of a date to ensure it matches the format
    YYYY/MM/DD. This function makes uses the of strptime function from
    the python datetime library. The strptime function tries to convert
    a date string to a datetime object, but requires a format. If the date
    string does not match the provided format, the function throws an error.
    Therefore, we can catch the error and use that to denote that the date
    is in incorrect format.  
    
    Returns true if the data is in the correct format or 
    false if it is not in the correct format.

    Parameters
    ----------
    date: str
        A date string to be checked.
    """

    # Wraps the strptime function in a try-catch block. The strptime
    # function throws an error if the date string is not in the expected
    # format. Therefore, if we catch an error, it means that the date
    # is not in the correct format.
    try:
        # Tries to run the strptime function using the format %Y/%m/%d which
        # is equivalent to YYYY/MM/DD
        dt.strptime(date, "%Y/%m/%d")
        
        # If this function completes successfully, return true.
        return True
    
    except ValueError:
        # If an error is caught, return false
        return False


def checkAccession(a):
    """Checks an record's accession number to see whether it contains
    a 'patent' sequence prefix. Patent sequences are things like vaccine or laboratory
    strains which may bias a tree or probe set (or whatever you want to use the output of this
    script for) toward a sequence that does not exist in nature. Returns true is the accession number does
    not contain one of these prefixes, and returns false if a patent sequence prefix is found.

    Parameter
    ---------
    a: str
        An accession number to check.
    """

    # Creates a boolean value to denote whether the accession number is valid.
    valid_accession = True

    # Uses regex to parse the accession number and capture the prefix
    r = re.search(accession_regex_pattern, a)

    # Checks whether a match was found. In testing, some cases have arisen
    # where protein sequences have been included (which do not have
    # the same accession numbers as nucleotide sequences) 
    if r:
                
        # If there was a match, the script then checks whether the captured prefix 
        # is found within the list of patent sequence prefixes.
        if r.group(1) in patent_sequence_prefixes:
            # If yes, the accession number is not valid and should
            # not be downloaded.
            valid_accession = False
    else:
        # If no match was found, meaning that the accession number did not match
        # a typical Nucleotide accession number, the accession number is not valid 
        # and should not be downloaded.
        valid_accession = False

    # Return whether the accession number is valid.
    return valid_accession


def checkCoverage(sequence, cov):
    """Checks whether the sequence is defined in NCBI and whether the sequence coverage
    (number of ambiguous bases (N's) divided by the total sequence length) and ensures it
    meets a minimum threshold. Returns true is the sequence passes and false is not.

    As a note, the Bio.Seq object's 'defined' method is an undocumented method that I found
    while searching the Biopython GitHub Issues (https://github.com/biopython/biopython/issues/3667). 
    It provides a convenient way of checking whether the record contains a genome sequence or not.

    Parameters
    ----------
    sequence: str
        A nucleic acid sequence to be checked

    cov: float
        The minimum percent coverage required for a sequence
        to pass
    """

    # Creates a boolean value to denote whether a sequence passed
    # the filters.
    valid_sequence = True

    # First, checks whether the sequence is defined. 
    if sequence.defined:
        # If so, check the coverage

        # Counts the number of N characters in a sequence using the 
        # count() function built into python.
        Ns = sequence.count("N")
    
        # Grabs the length of the sequence.
        seqLen = len(sequence)

        # Calculates percent coverage by taking the number of Ns divided by the 
        # total length of the sequence and subtracting from 1.
        coverage = (1 - (Ns / seqLen))

        # Checks whether the sequence coverage is below the coverage threshold
        if coverage < cov:
            # If so, then the sequence has failed the coverage filter and
            # the boolean is set accordingly.
            valid_sequence = False
    else:
        valid_sequence = False

    # Return whether the sequence coverage is valid.
    return valid_sequence

def checkSegment(recSegment, segLengths, recLen, lenVar):
    """Checks whether a sequence of a particular segment meets the 
    expected length range. The length range is calculated by multiplying 
    the average length of the segment (based on RefSeq) by a variance factor 
    and adding/subtracting this from the average to create a length range.
    Returns true is the segment passes and false if not.

    Parameters
    ----------
    recSegment: string
        The segment for a given record
    
    segLengths: dict
        A dictionary mapping segments to their average lengths
    
    recLen: int
        The length of a given record

    lenVar: float
        A percent variance to allow in the length of a sequence.
    """

    # Creates a boolean value to denote whether a segment passes 
    passing_segment = False

    # First check whether the sequence's segment is found within
    # the dictionary containing the average length of genome segments.
    # If the segment is not found, then it does not match the segments of
    # the RefSeq, and should be filtered out.
    if recSegment in segLengths.keys():

        # Next, calculate the minimum and maximum segment
        # lengths by multiplying the average segment length
        # by the variance factor and adding/subtracting this
        # from the average length.
        avg = segLengths[recSegment]
        minLen = avg - (avg * lenVar)
        maxLen = avg + (avg * lenVar)

        # Check whether the segment length falls within that range.
        if recLen >= minLen and recLen <= maxLen:
            # If so, then the sequence passes.
            passing_segment = True

    # Return whether the sequence passes
    return passing_segment

def writeErrorSummary(outPref, date, org, errorMsg):
    """A function that writes an error summary file in
    a set format, but allows various different errors to be
    passed to handle multiple error cases.

    Parameters
    ----------
    outPref: string
        A prefix to include in the file name in addition to 
        "-download-summary.txt".

    date: date
        The date that NCBI was accessed to include in the summary
        for reference by the user.

    org: string
        The organism supplied by the user that was queried.

    errorMsg: string
        The error message to be written to the file.
    """
    #Create a file stream to write the summary file.
    summaryOut = open("{0}-download-summary.txt".format(outPref), "w+")

    # Write data to the summary file including the date NCBI was accessed,
    # the organism, the genome length filters, the length variance,
    # the minimum coverage required, the header annotation parameters,
    # the submission date filter parameter, the query, the amount of sequences found,
    # the amount of sequences downloaded, and the amount of sequences filtered.
    summaryOut.write("Date Accessed: {0}\n".format(date))
    summaryOut.write("Organism: {0}\n".format(org))
    summaryOut.write("\n")
    summaryOut.write("Error Message: {0}\n".format(errorMsg))
    summaryOut.close()

def main():

    # Creates an argument parser and defines the script arguments
    parser = ap.ArgumentParser()

    parser.add_argument('-o', '--organism', required = True, type=str, \
        help='[Required] - Organism to be searched', \
        action = 'store', dest = 'org')
    parser.add_argument('-e', '--email', required = True, type=str, \
        help='[Required] - email address to be used to query NCBI', \
        action = 'store', dest='email')
    parser.add_argument('--output', required=True, type=str, \
        help='[Required] - output file prefix', \
        action='store', dest='outpref')
    parser.add_argument("--annotateHeader", required=False, \
        help = "Modifies the sequence id of downloaded genomes to contain the taxid and segment/chromosome before the accession number (Example: taxid|segment|accession)", \
        action='store_true', dest='annotHeader')
    parser.add_argument("--allowedLengthVariance", required = False, type=float, \
        help = 'The script will determine the organism\'s genome size by looking at the ref seq. This parameter sets the minimum and maximum sequence length for querying Genbank [Default: 0.05]', \
        action = 'store', dest='lenVar')
    parser.add_argument("--minCov", required=False, type=float, \
        help = "The minimum coverage required for a sequence to pass filters and be downloaded. [Default = 0.95]", \
        action='store', dest='minCov')
    parser.add_argument("--minLen", required=False, type=int, \
        help = '[Requires --maxLen parameter] - If supplied, the script will use the provided minimum and maximum length parameters to query NCBI. (NOTE - this will not consider any segments/chromosomes that fall outside of the length range provided) [Default = OFF]', \
        action = 'store', dest='minLen')
    parser.add_argument("--maxLen", required=False, type=int, \
        help = '[Requires --minLen parameter] - If supplied, the script will use the provided minimum and maximum length parameters to query NCBI. (NOTE - this will not consider any segments/chromosomes that fall outside of the length range provided) [Default = OFF]', \
        action='store', dest='maxLen')
    parser.add_argument("--separateSegments", required=False, \
        help = "If supplied, the script will output a separate file containing sequences for each segment/chromosome of the organism", \
        action='store_true', dest='sepOut')
    parser.add_argument("--minDate", required=False, type=str, \
        help = "[Requires --maxDate parameter] Only download sequences submitted after this date. (Must be in YYYY/MM/DD format) [Default = OFF]", \
        action='store', dest='minDate')
    parser.add_argument("--maxDate", required=False, type=str, \
        help="[Requires --minDate parameter] Only download sequences submitted before this date. (Must be in YYYY/MM/DD format) [Default = OFF]", \
        action='store', dest='maxDate')

    args = parser.parse_args()

    # Supplies the provided email to the NCBI entrez package.
    Entrez.email = args.email

    # Sets the minimum sequence coverage at 95% by default
    minCov = 0.95
    # If the minimum coverage option was supplied by the user, set
    # this value as the minimum sequence coverage instead.
    if args.minCov:
        minCov = args.minCov

    # Grabs the current date of accession (for final summary report at the end of the analysis).
    currentDate = dt.now().strftime("%Y/%m/%d")

    # We can first determine whether the provided organism has multiple segments/chromosomes
    # in its genome by querying the NCBI genome database.
    numSegments = 0
    genomeQuery = "({0}[Orgn])".format(args.org)
    genomeSearch = Entrez.esearch(db="genome", term=genomeQuery, retmax=retmax)
    # Parse the esearch result.
    genomeRecs = Entrez.read(genomeSearch)

    if len(genomeRecs['IdList']) == 0:
        error = "No NCBI Genome entry found for organism {0}".format(args.org)
        writeErrorSummary(args.outpref, currentDate, args.org, error)
        sys.exit(error)

    # Loop over each record in the esearch output (Should only be one.)
    for id in genomeRecs["IdList"]:
        # Use the esummary function to grab a summary for that record and
        # parse the result.
        idHandle = Entrez.esummary(db='genome', id=id, retmode='txt')
        idData = Entrez.read(idHandle)[0]

        # Grab the number of segments from the summary data.
        numSegments = int(idData['Number_of_Chromosomes'])

        if numSegments == 0:
            error = "No number of Chromosomes/Segments listed in the Genome entry for organism {0}".format(args.org)
            writeErrorSummary(args.outpref, currentDate, args.org, error)
            sys.exit(error)


    # Print a message to the console letting the user know how many
    # segments were identified. 
    print("{0} found to have {1} segment(s)/chromosome(s)\n".format(args.org, numSegments))

    ## Handles the date range options
    # We will create a string to be appended to a later query filtering by
    # date if a minimum/maximum date were supplied. For now, create an empty
    # string.
    queryDateString = ""
    
    # Check which options the user supplied
    if args.minDate and not args.maxDate:
        # If the user supplied a minimum date without a maximum date, notify the user and exit
        sys.exit("ERROR: --minDate option supplied without the --maxDate option. Please supply both options.")
    elif args.maxDate and not args.minDate:
        # If the user supplied a maximum date without a minimum date, notify the user and exit
        sys.exit("ERROR: --maxDate option supplied without the --minDate option. Please supply both options.")
    elif args.minDate and args.maxDate:
        # If the user supplied both a minimum and maximum date, then we can continue.
        
        # Checks whether the minimum date was provided in the correct format
        if not checkDateFormat(args.minDate):
            # If not, notify the user and exit
            sys.exit("ERROR: --minDate {0} is not in the correct format (YYYY/MM/DD). Please supply it in the correct format".format(args.minDate))

        # Checks whether the maximum date was provided in the correct format.
        if not checkDateFormat(args.maxDate):
            # If not, notify the user and exit.
            sys.exit("ERROR: --maxDate {0} is not in the correct format (YYYY/MM/DD). Please supply it in the correct format.".format(args.maxDate))

        # If both dates were provided in the correct format, we can construct the 
        # filter string to appended to the NCBI query.
        queryDateString = " AND (\"{0}\"[PDAT]:\"{1}\"[PDAT]))".format(args.minDate, args.maxDate)

    ## Handles the minimum and maximum sequence length options.
    # Set both to 0 by default.
    minLen = 0
    maxLen = 0

    # The length variance should be set at 5% by default
    lenVar = 0.05
    # Checks which options the user supplied. 
    if args.minLen and not args.maxLen:
        # If the user supplied a minimum length without a maximum length, notify the user and exit
        sys.exit("ERROR: --minLen option supplied without the --maxLen option. Please supply both options.")
    elif not args.minLen and args.maxLen:
        # If the user supplied a maximum length without a minimum length, notify the user and exit
        sys.exit("ERROR: --maxLen option supplied without the --minLen option. Please supply both options.")
    elif args.minLen and args.maxLen:
        # If the user supplied both a minimum and maximum length, then we can continue.

        # If the minimum and maximum length option were supplied, then the length variance option
        # will not have any affect on the analysis. Thus, we check whether the user has set it.
        if args.lenVar:
            # If length variance has been set, notify the user and exit (that way they know it has no effect on the analysis.)
            print("NOTE: --lenVar option was supplied, but a defined length range was also supplied. This parameter will be ignored.")

        # Set the minimum and maximum length variables to the values provided by the user.
        minLen = args.minLen
        maxLen = args.maxLen
    else:
        # If the user did not provide minimum and maximum length values, then we need to calculate
        # this range. This is done by querying NCBI for the organism's RefSeq sequence and taking
        # the average sequence length. If no RefSeq exists for that organism, then the keyword
        # 'complete genome' is added to the query, and the lengths of the return sequences are averaged.

        # If the user specified a length variance value using the option,
        # set that value 
        if args.lenVar:
            lenVar = args.lenVar

        # Query NCBI's nucleotide database for the organism's RefSeq genome sequence(s).
        refSeqQuery = "({0}[Orgn]) AND refseq[filter]".format(args.org)
        refSeqSearch = searchNucleotide(refSeqQuery)
        seqTypeFound = "RefSeq"

        # Check whether any record IDs were returned.
        if len(refSeqSearch["IdList"]) == 0:
            # If not, the script will use different title qualifiers to return
            # sequences which likely represent a subset of full genome sequences.
            # The list of qualifiers is globally defined for easy modification.
            print("No RefSeq sequences are defined for organism {0}...\n".format(args.org))
            seqTypeFound = "Complete Genome"
            for qual in title_search_quantifiers:
                print("Trying sequences containing '{0}' in the title...".format(qual))
                query = "({0}[Orgn]) AND {1}[title]".format(args.org, qual)
                print(query)
                refSeqSearch = searchNucleotide(query)

                if len(refSeqSearch["IdList"]) != 0:
                    break
                else:
                    print("No sequences found with titles matching that criteria...\n")

            if len(refSeqSearch["IdList"]) == 0:
                error = "No RefSeq or complete genome sequences found for organism {0}".format(args.org)
                writeErrorSummary(args.outpref, currentDate, args.org, error)
                sys.exit(error)

            #print("No RefSeq exists for {0}... Searching for sequences containing 'complete genome' in the title\n".format(args.org))
            #completeGenomeQuery = "({0}[Orgn]) AND complete genome[title]".format(args.org)
            #refSeqSearch = searchNucleotide(completeGenomeQuery)


        # Create empty variables to store the total length of all sequences returned and 
        # the number of records. Because there could be multiple segments/chromosomes, a dictionary
        # is used to store these values. If the organism contains only a single segment/chromosome,
        # 0 will be used as the key in the dictionary.
        #totalLength = {}
        segLens = {}
        numRecords = {}

        # Loop from 0 to the number of IDs returned from the RefSeq query increasing by 10,000
        # at each step. Entrez is optimized for batch searches of 10,000 records at a time. Thus, this will
        # optimize the script by limiting the number of times we ping NCBI.
        for i in range(0, len(refSeqSearch['IdList']), 9999):

            # Fetch the records from the nucleotide database.
            idHandle = Entrez.efetch(db='nucleotide', id = ",".join(refSeqSearch["IdList"]), retstart=i, retmax=retmax, rettype='gb', retmode="txt")
            
            # Loop over the GenBank records returned from NCBI.
            for seqRecord in SeqIO.parse(idHandle, "genbank"):

                # Grab the length of the sequence
                # and the accession number.
                seqLength = len(seqRecord.seq)
                accession = seqRecord.id
            
                # Check the accession number of the sequence to ensure that
                # it is not a patent sequence (this could skew the average length as these
                # sequences may not be representative of the actual viral genome).
                if checkAccession(accession):

                    # If the accession number passes filtering, then we can 
                    # continue processing.

                    # Grab the source feature and description from the seqRecord
                    recSource = getSourceFeature(seqRecord)
                    recDesc = seqRecord.description

                    # Grab the records' segment using the source and 
                    # description.
                    segment = getRecSegment(recSource, recDesc)

                    # Checks whether multiple segments noted in the genome entry for the organism,
                    # but the segment could not be determined for the given RefSeq.
                    if segment == None and numSegments != 1:
                        # If this is the case, set the segment to 0
                        segment = 0

                        ## TO BE ADDED - using named segments instead of numbers

                        #print("No Segment Number Annotation found for Organism. First gene product found on segment will be used to name segment.")

                        #product = getFirstGeneProduct(seqRecord)

                        #if product == None:
                        #    error = "No Segment Number or gene product annotation found to differentiate viral segments for organism {0}".format(args.org)
                        #    writeErrorSummary(args.outpref, currentDate, args.org, error)
                        #    sys.exit(error)
                        #else:
                        #    segment = product

                    # If no segment was found in the RefSeq entry, but the organism
                    # only contains 1 segment, then we can simply set the segment equal to 0.
                    elif segment == None and numSegments == 1:
                        segment = 0
            
                    # Checks whether the segment was found for an organism with 
                    # a segmented genome or whether the organism only has 1 segment.
                    if (numSegments > 1 and segment != 0) or (numSegments == 1):
                        # If so, we can add the length of the current RefSeq to a 
                        # dictionary to be used in calculating the average length of that organism's genome.

                        # If no data for the current segment has been added to the
                        # totalLength dictionary, create a new entry for it.
                        if segment not in segLens.keys():
                            segLens[segment] = []

                        # If no data for the current segment has been added to the
                        # numRecords dictionary, create a new entry for it.
                        if segment not in numRecords.keys():
                            numRecords[segment] = 0

                        # Increment the number of records for the segment by 1
                        # and add the length of this segment to the total length 
                        numRecords[segment] += 1
                        segLens[segment].append(seqLength)


        # One possibility for error at this point is if a RefSeq does not
        # exist for all of the possible segments of an organism. If so,
        # the script will exit and write an error summary file notifying hte user.
        if len(segLens.keys()) != numSegments:
            error = "ERROR: Number of segments returned from RefSeq does not match number reported in genome entry for organism {0}".format(args.org)
            writeErrorSummary(args.outpref, currentDate, args.org, error)
            sys.exit(error)

        # Create an empty dictionary to store the average lengths for
        # each segment/chromosome (again, if organism only has 1 segment/chromosome
        # then it will be under segment '0').
        avgLengths = {}

        # Create variables to store the minimum and maximum
        # average lengths. This is done because, if the organism 
        # has multiple segments/chromosomes, we want to capture all of the segments
        # in our query. Thus, the length range will need to contain all of the segments/chromosomes.
        # The average length of the smallest segment will be the minimum and the average length
        # of the largest segment will be maximum. If only a single segment/chromosome,
        # the values will be the same (as we only need to capture this one segment)
        minAvg = 0
        maxAvg = 0

        # Calculating the length range for the Nucleotide query
        # is handled differently depending on whether the organism contains more than 1 
        # segment or not.
        if len(segLens.keys()) > 1:
            # If the organism contains multiple segments, then the minimum genome length
            # should be calculated by taking the average length of the smallest segment
            # minus the length variance. The maximum genome length should be calculated as
            # by taking the average length of the largest segment plus the length variance.
            
            # Loop over each segment in the total length dictionary.
            for seg in sorted(segLens.keys()):

                # If no data for the current segment has been added to the avgLengths
                # dictionary, create a new entry for it.
                if seg not in avgLengths.keys():
                    avgLengths[seg] = 0
        
                # Calculate the average length by dividing the total length
                # by the number of records for the current segment.
                avgLengths[seg] = sum(segLens[seg]) / numRecords[seg]

                # If the minimum average has not been set yet (is equal to 0) or
                # the average length is less than the current minimum. Set
                # the minimum average to the current average length
                if minAvg == 0 or avgLengths[seg] < minAvg:
                    minAvg = avgLengths[seg]

                # If the maximum average has not been set yet (is equal to 0) or
                # the average length is greater than the current maximum. Set
                # the maximum average to the current average length
                if maxAvg == 0 or avgLengths[seg] > maxAvg:
                    maxAvg = avgLengths[seg]

                print("Segment {0}: {1} {2} Sequences\nAverage Length {3} bp".format(seg, numRecords[seg], seqTypeFound, avgLengths[seg]))
        
        else:
            # If the organism only contains 1 segment in the genome,
            # then the length range will be the smallest and largest RefSeq lengths 
            # +/- the length variance factor.

            # Grab only segment key from the dictionary (this is
            # added for the future in case any formal segment/gene names will be used)
            seg = list(segLens.keys())[0]

            # Grab the minimum and maximum RefSeq genome lengths for
            # the organism from the dictionary.
            minAvg = min(segLens[seg])
            maxAvg = max(segLens[seg])

            print("{0} {1} Sequences\nLargest Sequence: {2} bp\nSmallest Sequence: {3} bp".format(numRecords[seg], seqTypeFound, maxAvg, minAvg))
            


            # Print the number of RefSeq Sequences and average length to the console
            # for the user.
            #if seg != 0:
                # If the segment does not equal 0 (the organism has segments/chromosomes), our message will
                # reflect this by also printing the segment
            #    print("Segment {0}: {1} RefSeq Sequences\nAverage Length {2}".format(seg, numRecords[seg], avgLengths[seg]))
            #else:
                # If the segment is 0, then it means this organism does not have any segments/chromosomes,
                # and we only need to print the number of sequences and average length.
            #    print("{0} RefSeq Sequences\nAverage Length {1}".format(numRecords[seg], avgLengths[seg]))
    
        # Print a blank line to separate text in the console.
        print("")
    
        # Calculates the minimum length of query by adjusting the minimum
        # average by the length variance factor.  
        minLen = minAvg - (minAvg * lenVar)

        # Calculates the maximum length of query by adjusting the maximum
        # average by the length variance factor.  
        maxLen = maxAvg + (maxAvg * lenVar)

    # Create the search term containing the organism, length range, and date filter (will be blank unless
    # supplied by the user).
    searchTerm = "({0}[Orgn]) AND {1}:{2}[SLEN]{3}".format(args.org, minLen, maxLen, queryDateString)

    # Print the query to console for the user.
    print("NCBI Query: {0}\n".format(searchTerm))
            
    # Query NCBI nucleotide using the history option (optimization for large queries),
    # parse the output, and close the record handle.
    recordsHandle = Entrez.esearch(db="nucleotide", term=searchTerm, retmax=retmax, usehistory="y")
    searchResults = Entrez.read(recordsHandle)
    recordsHandle.close()

    # Grab the list of Ids returned by the query and count the
    # number returned.
    idsFound = searchResults['IdList']
    total_seqs = len(idsFound)

    # Query the ids returned from teh initial esearch using the
    # EPOST function.
    ePostHandle = Entrez.epost("nucleotide", id = ",".join(idsFound))
    ePostResult = Entrez.read(ePostHandle)

    # Grab the WebEnv and Query key for the history feature (optimization for large queries)
    webEnv = ePostResult["WebEnv"]#searchResults["WebEnv"]
    query_key = ePostResult["QueryKey"]#searchResults["QueryKey"]

    # Print a message letting the user know how many records were found
    print("Query returned {0} Sequences\n\n".format(total_seqs))

    # Print the message letting the user know what filtering will be applied.
    if numSegments > 1 and (not args.minLen and not args.maxLen):
        # If there are multiple segments and the user did not supply a set length range, we will perform
        # filtering based on the accession length, coverage percentage, and length of the segment.
        print("Filtering any 'patent' sequence by accession prefix, sequences with less than {0}% coverage, or segment length outside the segment specific range\n".format(str(minCov * 100)))
    else:
        # Otherwise, we will simply filter based on accession length and coverage percentage
        print("Filtering any 'patent' sequence by accession prefix or sequences with less than {0}% coverage\n".format(str(minCov * 100)))

    # Create an empty dictionary to contain output streams for each segment. If
    # the user supplied the --separateSegment option and the genome has multiple
    # segments/chromosomes, this will be used. If not, then a single output stream will
    # be created at key value 0
    seqOut = {}

    # Create filestreams for metadata from passing and filtered sequences.
    metadataOut = open("{0}-metadata.tsv".format(args.outpref), "w+")
    filteredOut = open("{0}-filtered-metadata.tsv".format(args.outpref), "w+")

    # Write headers to the metadata files.
    if numSegments > 1:
        # If organism has multiple segments/chromosomes, then include a
        # column for the segment in the file.
        metadataOut.write("Accession Number\tLength\tTaxonomy ID\tSegment\tDescription\tIsolate\tContinent\tCountry\tRegion\tCollection Date\tIsolate Source\tAuthors\tTitle\tJournal\n")
        filteredOut.write("Accession Number\tLength\tTaxonomy ID\tSegment\tDescription\tIsolate\tContinent\tCountry\tRegion\tCollection Date\tIsolate Source\tAuthors\tTitle\tJournal\n")
    else:
        # If the organims has a single segment/chromosome, do not 
        # include a column for segment (as it does not exist)
        metadataOut.write("Accession Number\tLength\tTaxonomy ID\tDescription\tIsolate\tContinent\tCountry\tRegion\tCollection Date\tIsolate Source\tAuthors\tTitle\tJournal\n")
        filteredOut.write("Accession Number\tLength\tTaxonomy ID\tDescription\tIsolate\tContinent\tCountry\tRegion\tCollection Date\tIsolate Source\tAuthors\tTitle\tJournal\n")


    # Create variables to store the number of sequences downloaded vs
    # filtered.
    downloaded_seqs = 0
    filtered_seqs = 0
 
    count = 0

    #print(idsFound[9999])
    #exit()

    # Loop from 0 to the number of IDs returned from the RefSeq query increasing by 10,000
    # at each step. Entrez is optimized for batch searches of 10,000 records at a time. Thus, this will
    # optimize the script by limiting the number of times we ping NCBI.
    for i in range(0, len(idsFound), 9999):

        # Subset the list of Ids to pull using the current iterator
        idsToGrab = idsFound[i:i+10000]

        # Use NCBI's efetch function with the ePOST webenv/query key to 
        # Grab the Genbank records for these Ids
        idHandle = Entrez.efetch(db='nucleotide', id = ",".join(idsToGrab), retmax=10000, rettype='gb', retmode="txt", webenv=webEnv, query_key=query_key)

        # Parse the returned handle into SeqRecords
        seqRecords = SeqIO.parse(idHandle, "genbank")

        # Including this in a try-catch block to try and catch "incomplete read"
        # exceptions which occur when communication with NCBI was interrupted due
        # to connection.
        #
        # Have yet to figure out a concrete solution besides running the script on
        # ethernet.
        try:

            # Loop over each ID in the list list retrieved.
            for seqRecord in seqRecords:
                count += 1

                # Grab the accession number and description
                # from the record
                accession = seqRecord.id
                print("Processing Record: {0}, Number {1} of {2}, Sequence Defined = {3}".format(accession, count, len(idsFound), seqRecord.seq.defined))
                desc = seqRecord.description

                # Get the source feature from the record.
                recSource = getSourceFeature(seqRecord)

                # Grab the first reference from the record's annotations
                recRef = seqRecord.annotations['references'][0]

                # Parse the source feature and reference to gather metadata 
                # values (if provided)
                taxID, isolate, continent, country, region, iso_source, col_date, authors, title, journal = getReportData(recSource, recRef)

                # Grab the length of the the sequence.
                recLen = len(seqRecord.seq)
                if seqRecord.seq.defined == False:
                    recLen = "Sequence Undefined in NCBI"

                # Create an empty list to store the metadata to be written for each sequence
                outMetadata = []
                
                # Also create a boolean variable to denote whether a sequence has passed filters or not. 
                seqPass = False

                # Check the number of segments/chromosomes that the organism should have.
                if numSegments > 1:

                    # If more than 1, then we will output the segment/chromosome 
                    # to the metadata.
                    
                    # Grab the segment for the record.
                    recSeg = getRecSegment(recSource, desc)

                    #if recSeg == None:
                    #    recSeg = getFirstGeneProduct(seqRecord)

                    # Create an empty variable to store the name of the output
                    # filestream (used for if the user specified to separate
                    # sequences by segment/chromosome). If the user did not specify
                    # this option, then this value will remain blank and all sequences will be written
                    # to the same file.
                    outStream = ""

                    # Check whether the user specified the separate output option. 
                    if args.sepOut:
                        # If yes, determine the output stream to write the sequence to.
                        if recSeg == None:
                            # If the segment is None, then the segment will be written to an
                            # unknown file
                            outStream = "unknown-segment-"
                        else:
                            # Create the output stream named segment-{segment value}-
                            # (this will get concatenated into the output file name)
                            outStream = "segment-{0}-".format(recSeg)

                    # If the record segment is none, change the value to unknown for writing to metadata.
                    if recSeg == None:
                        recSeg = "Unknown"

                    # Compile the sequence's metadata into a list
                    outCols = [accession, str(recLen), taxID, str(recSeg), desc, isolate, continent, country, region, col_date, iso_source, authors, title, journal]

                    ## Perform filtering

                    # Check whether the user provided a custom minimum and maximum length range
                    if args.minLen and args.maxLen:
                        # If yes, filter by checking the accession number and coverage (ignore
                        # per segment length filtering)
                        if checkAccession(accession) and checkCoverage(seqRecord.seq, minCov):
                            # If the sequence passes the accession and coverage checks, mark
                            # it as passing
                            seqPass = True
                    else:
                        # If no, filter by checking the accession number, coverage, and segment
                        # lengths.
                        if checkAccession(accession) and checkCoverage(seqRecord.seq, minCov) and checkSegment(recSeg, avgLengths, recLen, lenVar):
                            # If the sequence passes the checks, mark it as passing
                            seqPass = True

                    # If the sequence passed filters and the user supplied the annotate header
                    # option, modify the sequence header to include the taxID, segment, and 
                    # accession number.
                    if args.annotHeader and seqPass == True:
                                seqRecord.id = "{0}|{1}|{2}".format(taxID, accession, "segment_" + recSeg)
                else:

                    # If a single segment/chromosome, no segment/chromosome
                    # will be written to metadata.

                    # Create an empty variable to store the name of the output
                    # filestream (used for if the user specified to separate
                    # sequences by segment/chromosome). If the user did not specify
                    # this option, then this value will remain blank and all sequences will be written
                    # to the same file.
                    outStream = ""

                    # Compile the sequence's metadata into a list
                    outCols = [accession, str(recLen), taxID, desc, isolate, continent, country, region, col_date, iso_source, authors, title, journal]

                    ## Perform filtering

                    # Filter by checking the accession number and coverage of
                    # the sequence.
                    if checkAccession(accession) and checkCoverage(seqRecord.seq, minCov):
                        # If the sequence passes the accession number and coverage checks, 
                        # mark it as passing.
                        seqPass = True

                        # If the user supplied the annotate header option, modify
                        # the sequence header to include the taxID, "complete genome", and
                        # accession number.
                        if args.annotHeader:
                            seqRecord.id = "{0}|{1}|{2}".format(taxID, accession, "complete_genome")
                
                # Check whether the sequence passed filtering
                if seqPass:
                    # If yes, we will write the sequence to the filestream determined earlier.

                    # If the file stream is not part of the seqOut dictionary (i.e. is has not 
                    # been created yet), open a new output fasta file.
                    if outStream not in seqOut.keys():
                        seqOut[outStream] = open("{0}-{1}sequences.fasta".format(args.outpref, outStream), "w+")

                    # Write the sequence to the fasta output stream
                    SeqIO.write(seqRecord, seqOut[outStream], 'fasta')

                    # Write the metadata to the passing metadata file
                    metadataOut.write("\t".join(outCols) + "\n")

                    # Increase the number of downloaded sequences by 1
                    downloaded_seqs += 1 
                else:
                    # If no, write the metadata to the filtered metadata file. 
                    filteredOut.write("\t".join(outCols) + "\n")

                    # Increase the number of filtered sequences by 1
                    filtered_seqs += 1

        # If an Incomplete Read error occurs, print a message to console letting the user know
        except (http.client.IncompleteRead) as e:
            print("Error reading batch")

        # Once all of the GenBank records have been parsed, close the efetch handle.
        idHandle.close()


    # Print messages to the console letting the user know how many sequences were
    # downloaded and filtered.
    print("Downloaded {0} sequences".format(downloaded_seqs))
    print("Filtered {0} sequences".format(filtered_seqs))

    
    # Loop over the output file streams and close them.
    for k in seqOut.keys():
        seqOut[k].close()

    # Close the metadata file streams.
    metadataOut.close()
    filteredOut.close()


    ## Create the summary file

    # Convert the length variance factor into a string with 
    # a percent symbol for the summary file.
    lenVarSummary = str((lenVar * 100)) + "%"
    # If the user supplied a minimum and maximum length,
    # then change the length variance summary value to
    # a message letting the user know that length variance
    # was disabled.
    if args.minLen and args.maxLen:
        lenVarSummary = "Length range provided by user - length variance disabled."

    # Create a variable to summarize the date filter parameters. The
    # filter is disabled by default.
    dateFilter = "Disabled"
    # If the user supplied and minimum and maximum date, update the
    # summary variable to include the date range in string format.
    if args.minDate and args.maxDate:
        dateFilter = "{0}-{1}".format(args.minDate, args.maxDate)

    # Create a variable to summarize the header annotation
    # parameter. The option is disabled by default.
    headerAnnotate = "Disabled"
    # If the user supplied the annotate header option, 
    # update the summary to note that this was enabled,
    if args.annotHeader:
        headerAnnotate = "Enabled"
    

    # Create a file stream to write the summary file.
    summaryOut = open("{0}-download-summary.txt".format(args.outpref), "w+")

    # Write data to the summary file including the date NCBI was accessed,
    # the organism, the genome length filters, the length variance,
    # the minimum coverage required, the header annotation parameters,
    # the submission date filter parameter, the query, the amount of sequences found,
    # the amount of sequences downloaded, and the amount of sequences filtered.
    summaryOut.write("Date Accessed: {0}\n".format(currentDate))
    summaryOut.write("Organism: {0}\n".format(args.org))
    summaryOut.write("Genome Length Range: {0}-{1}\n".format(minLen, maxLen))
    summaryOut.write("Allowed Genome Length Variance: {0}\n".format(lenVarSummary))
    summaryOut.write("Minimum Genome Coverage: {0}\n".format(str(minCov * 100) + "%"))
    summaryOut.write("Header Annotation: {0}\n".format(headerAnnotate))
    summaryOut.write("Submission Date Range: {0}\n".format(dateFilter))
    summaryOut.write("\n")
    summaryOut.write("Query: {0}\n".format(searchTerm))
    summaryOut.write("Total Sequences Returned: {0}\n".format(total_seqs))
    summaryOut.write("Sequences Filtered: {0}\n".format(filtered_seqs))
    summaryOut.write("Sequences Downloaded: {0}\n".format(downloaded_seqs))

    # Close the summary file.
    summaryOut.close()

if __name__ == "__main__":
    main()