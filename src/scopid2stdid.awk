#!/bin/awk -f
#
# Usage ./scopid2stdid.awk scopid-list-file
#
# if the scopid start with 'g', meaning a genetic domain sequence, output the
# pdbid instead
#
# created 2007-08-24, updated 2010-01-18, Nanjiang
#
# Copyright© Nanjiang Shu 
# Department of Materials and Environmental Chemistry, Stockholm University,
# Sweden
# Email: nanjiang.shu@mmk.su.se

function usage(e1)
{
    print "Usage:  scopid2stdid.awk [-p] [files...]"
    print " -p|-print-scopid  : print scopid as well"
    print ""
    print "Created 2007-08-24, updated 2010-01-18, Nanjiang Shu"
}
BEGIN {
# How many lines
    isPrintSCOPID = 0;
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
            isPrintSCOPID=1;
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
    pdbid = toupper(substr($1,2,4))
    firstchar=substr($1,1,1)
    if (firstchar == "d" )
    {
        chainID = substr($1,6,1);

    }
    else if(firstchar == "e")
    {
        chainID = substr($1,8,1);
    }
    chainID=toupper(chainID);
    if (chainID == "_")
    { # for the SCOP version ealier than 1.73, the chainID of the domain sequences from some PDB files with the " " chainID is set to "_", 2010-01-18 
        sysstring=sprintf("getpdbchainidlist %s",pdbid );
        sysstring | getline chainIDList;
        chainID=substr(chainIDList,1,1);
    }

    if (isPrintSCOPID == 1)
    {
        printf("%s ", $1);
    }
    if (chainID == "_" ||firstchar == "g")
    {
        printf("%s\n", pdbid);
    }
    else
    {
        printf("%s%s\n", pdbid,chainID);
    }
}
END{
    }
