#!/bin/sh

# Snippet for loading file from stream or from path as zipped or not zipped file

# USAGE: 
# toolName file

#############################################
###### START line80 ################################
# read the mutation file followed by the stream data (pileup of tumor/tumor+PON)
mawk '
BEGIN { 
    # set the file var or set to stream
    file = "'${1-"stream"}'";
    # allow for writing the word "stream" as arg
    file = (tolower(file) == "stream") ? "stream" : file;

    # setting read command depending on file extension
    if (file != "stream") {
        if (file ~ /.gz$/) {
            readCMD = "cat " file " 2>/dev/null | gunzip ";
        } else {
            readCMD = "cat " file " 2>/dev/null ";
        } 
        # open file getline stream (Gstream) and skip first line
        # getline == 0 if file not found --> PON > "none
        if ((readCMD | getline) == 0) {
            printf("<toolName> File %s not found!\n", file) > "/dev/stderr";
            exit;
        } else {                   
            printf("<toolName> Reading from file %s.\n", file) > "/dev/stderr";
            getFile=1; # stateVariable to activate getline for file
            # not needed here but good pattern for combination of stream and fileRead
        }
        # if there is no stream data everything has to happen in BEGIN {}
        if (getFile == 1) {
            while ((readCMD | getline) > 0) {
                print($0);
            }
            exit;
        }
    } else {
            printf("<toolName> Reading from stream\n") > "/dev/stderr";
        }
}
# read stream data
{ # if readFile is set
    print($0);
}'