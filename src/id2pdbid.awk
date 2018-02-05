#!/bin/awk -f
#
# Usage  id2pdbid.awk  files
#
# extract the pdbid from the normal ids
# ids can be of 
# 1. five or four letter PDB chain id, e.g. 1BMFA 1AOL
# 2. scopid
#
# created 2010-01-17, updated 2010-05-27, Nanjiang
#
# Copyright© Nanjiang Shu 
# Department of Materials and Environmental Chemistry, Stockholm University,
# Sweden
# Email: nanjiang.shu@mmk.su.se

function usage(e1)
{
    print "Usage:  id2pdbid.awk [-p] [files...]"
    print "   ids can be "
    print "   1. five or four letter PDB chain id, e.g. 1BMFA 1AOL"
    print "   2. scopid, e.g. d1a12a_ g1o7d.3 "
    print " -p|-print-id  : print the original input id as well"
    print ""
    print "Created 2010-01-17, updated 2010-05-27, Nanjiang Shu"
    print ""
    print "Examples:"
    print ""
    print "   id2pdbid.awk idlist.txt"

}
BEGIN {
# How many lines
    isPrintID = 0;
    isExit = 0;
    i = 1
    while (i < ARGC) 
    {
        if (ARGV[i] == "-h" || ARGV[i] == "--help")
        {
            isExit = 1;
            usage();
            exit 1;
        }
        else if (ARGV[i] == "-p" || ARGV[i] == "--print-scopid")
        {
            isPrintID=1;
            delete ARGV[i];
            i ++;
        }
        else
        {
            i++;
        }
    }
}
{
    lengthID = length($1);
    firstchar = substr($1,1,1);
    pdbid="UNK"
    if (lengthID >=4 && lengthID <=5 )
    {
        pdbid = tolower(substr($1,1,4));
    }
    else if (lengthID >=7 && lengthID <=8 )
    {
        if(firstchar == "d" || firstchar == "e" || firstchar == "g")
        {
            pdbid = substr($1,2,4)
        }
    }
    if (pdbid=="UNK")
    {
        print "Error! Unrecognized id ", $1, " at line", NR  >"/dev/stderr"
    }
    if (isPrintSCOPID == 1)
    {
        print $1, pdbid
    }
    else 
    {
        print pdbid
    }
}
END{
    }
