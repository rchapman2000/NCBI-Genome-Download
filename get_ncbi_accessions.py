import sys
import os
import re
import json
import http
import argparse as ap

from Bio import Entrez, SeqIO 
import pycountry_convert as pc
from datetime import datetime as dt

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

def main():

    # Creates an argument parser and defines the script arguments
    parser = ap.ArgumentParser()

    parser.add_argument("-i", "--input", required = True, type=str, \
        help = '[Required] - a file containing a list of accession numbers to download from NCBI', 
        action = 'store', dest = 'inFile')
    parser.add_argument('-e', '--email', required = True, type=str, \
        help='[Required] - email address to be used to query NCBI', \
        action = 'store', dest='email')
    parser.add_argument('-o', '--outpref', required=True, type=str, \
        help='[Required] - output file prefix', \
        action='store', dest='outPref')
    
    args = parser.parse_args()

    # Supplies the provided email to the NCBI entrez package.
    Entrez.email = args.email

    # Handles the input parameter
    inFile = ''
    # Checks whether the file supplied by the user exists.
    if not os.path.exists(args.inFile):
        # If not, notify the user that the file does not exist and exit.
        sys.exit("ERROR: File {0} does not exist. Please supply an existing file.".format(args.inFile))
    else:
        # If the file does exist, open the file with read access
        inFile = open(args.inFile, "r")
    
    # Create an empty list to store accession numbers to search.
    accessions = []

    # Loop over each line in the input file and add
    # the accession number to the list.
    for line in inFile:
        accessions.append(line.strip("\n"))

    # Close the input file stream
    inFile.close()

    # Create an output file stream for sequence metadata and write a header.
    metadataOut = open("{0}-metadata.tsv".format(args.outPref), "w+")
    metadataOut.write("Accession Number\tLength\tTaxonomy ID\tDescription\tIsolate\tContinent\tCountry\tRegion\tCollection Date\tIsolate Source\tAuthors\tTitle\tJournal\n")

    # Create an output file stream for individual sequences.
    seqOut = open("{0}-sequences.fasta".format(args.outPref), "w+")

    # Loop over each accession number from the input file. A
    # while loop was used here to account for cases when a connection
    # to NCBI fails. This has been observed to occur when connecting via
    # a wifi connection. A while loop best handles this because we can leave the accession number
    # in the list if the connection fails, but can simply remove the accession number 
    # if the connection was successful.
    while len(accessions) != 0:

        # Grab the accession number from the start of the list.
        acc = accessions[0]

        # Print the accession number to the console.
        print(acc)

        # Wraps the NCBI connection code in a try-except block. If the connection to NCBI
        # is interrupted or dropped, an exception is thrown.
        try:
            # Search the Nucleotide database for the accession number
            # and parse the result.
            recHandle = Entrez.esearch(db="nucleotide", term=acc, retmax=1)
            searchResults = Entrez.read(recHandle)

            # Checks whether the search returned any results.
            if len(searchResults['IdList']) > 0:
                
                # If the search did return results, grab 
                # the first ID from the list (the top hit).
                recID = searchResults['IdList'][0]

                # Use the Entrez Efetch functionality to grab the Genbank
                # record for this ID.
                idHandle = Entrez.efetch(db='nucleotide', id = recID, rettype='gb', retmode="txt")

                # Parse the Genbank records returned using Biopython's SeqIO module
                seqRecords = SeqIO.parse(idHandle, "genbank")

                # Loop over the Genbank Records returned (there
                # should only be one)
                for seqRec in seqRecords:
                    # Grab the accession number and description
                    # from the record
                    accession = seqRec.id
                    desc = seqRec.description

                    # Get the source feature from the record.
                    recSource = getSourceFeature(seqRec)

                    # Grab the first reference from the record's annotations
                    recRef = seqRec.annotations['references'][0]

                    # Parse the source feature and reference to gather metadata 
                    # values (if provided)
                    taxID, isolate, continent, country, region, iso_source, col_date, authors, title, journal = getReportData(recSource, recRef)

                    # Grab the length of the the sequence.
                    recLen = len(seqRec.seq)

                    outCols = [accession, str(recLen), taxID, desc, isolate, continent, country, region, col_date, iso_source, authors, title, journal]

                    # Write the sequence to the fasta output stream
                    SeqIO.write(seqRec, seqOut, 'fasta')

                    # Write the metadata to the passing metadata file
                    metadataOut.write("\t".join(outCols) + "\n")

                    # Notify the user that the record was downloaded.
                    print("Record Downloaded\n")
            else:
                # If no record was returned for the accession number
                # provided, then the accession number must not exist. Notify
                # the user and exit.
                print("No NCBI record found for accession {0}\n".format(acc))

            # Remove the accession number from the list,
            # as it has been successfully processed
            accessions.pop(0)

        # If an Incomplete Read error occurs, print a message to console letting the user know
        except (http.client.IncompleteRead) as e:
            print("Connection to NCBI interrupted while querying for Accession {0}. Retrying...\n".format(acc))

    # Close output file streams
    metadataOut.close()
    seqOut.close()

if __name__ == "__main__":
    main()