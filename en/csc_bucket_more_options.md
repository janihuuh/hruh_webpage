Added March 1 2024 by Johannes

# 1. Get access to HRUH bucket at CSC Allas

To get started with the CSC buckets, please create a CSC account. https://my.csc.fi/login . Then ask Francesca to add you to the HRUH CSC project. 

The data files are encrypted, so you will also need the private key and the password to decrypt the data. You can get these from Francesca as well. 

Skip these steps if you already have these things.

# 2. Download data from a CSC bucket

To download data from a CSC bucket, you have at least these three options:

1.	Command line approach, which requires more IT expertise. Better for large files (hundreds of gigabytes). https://janihuuh.github.io/hruh_webpage/en/#!csc_bucket.md
2.	Using Cyberduck or similar tool to connect to CSC Allas, and then you can select the files that you want from the bucket. https://docs.csc.fi/data/Allas/using_allas/cyberduck/
3.	Go to https://sd-connect.csc.fi/ and download the files from there.  


# 3. Decrypt files

To decrypt the files, you can use this Graphical User Interface version of crypt4gh https://github.com/CSCfi/crypt4gh-gui (easier for a non-IT person) or the command line interface https://github.com/EGA-archive/crypt4gh
