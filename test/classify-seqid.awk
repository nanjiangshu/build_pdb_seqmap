#!/usr/bin/awk -f
#
# Usage   classify-seqid.awk files
#
# extract the ids with the specified idtype
# idtype can be of 
#  idtype = 0: five-letter pdb chain id
#  idtype = 1: scop domain id derived from single chain
#  idtype = 2: scop domain id (starts with "g") derived from multiple chains
#  idtype = 3: scop domain id (starts with "e")
#
# awk treats -h -t -m anyway as it's intrinsic arguments if it is argv[1]
# created 2010-01-17, updated 2015-02-20, Nanjiang
#
function usage(e1)
{
    print "Usage:   classify-seqid.awk  [files...] [-ldebug]"
    print " idtypes can be "
    print "   idtype = 0: five-letter pdb chain id"
    print "   idtype = 1: scop domain id derived from single chain"
    print "   idtype = 2: scop domain id (starts with "g") derived from multiple chains"
    print "   idtype = 3: scop domain id (starts with "e")"
    print " -ltype 0|1|2|3 : print the original input id as well"
    print " -lhelp         : print the help message and exit"
    print ""
    print "Created 2010-01-17, updated 2015-02-20, Nanjiang Shu"
    print ""
    print "Examples:"
    print ""
    print "  classify-seqid.awk -ltype 0 idlist.txt"
    print "  classify-seqid.awk -ltype 03 idlist.txt"
    print "# print a list of unrecognized ids"
    print "  classify-seqid.awk -ltype -1 idlist.txt"

}
BEGIN {
# How many lines
    idTypeList="";
    isDebug="false";
    i = 1;
    cntUsedARGC = 0;
    while (i < ARGC)
    {
        firstchar = substr(ARGV[i],1,1);
        if (firstchar == "-"){
            if (ARGV[i] == "-lhelp" )
            {
                isExit = 1;
                usage();
                delete ARGV[i]
                i ++;
                exit 1;
            }
            else if (ARGV[i] == "-ltype")
            {
                idTypeList=ARGV[i+1];
                cntUsedARGC +=2;
                delete ARGV[i]
                delete ARGV[i+1]
                i += 2;
            }
            else if ( ARGV[i] == "-ldebug")
            {
                isDebug="true";
                delete ARGV[i]
                cntUsedARGC ++;
                i++;
            }
            else
            {
                print "Bad argument", ARGV[i];
                exit 1;
            }
        }
        else
        {
            i++;
        }
    }
    sub("-1", "9", idTypeList);
}
{
    lengthID = length($1);
    firstchar = substr($1,1,1);
    idType = "9"; # unrecognized idType
    if (lengthID ==5 )
    {
        idType = "0";
    }
    else if (lengthID >=7 && lengthID <=8 )
    {
        if(firstchar == "d" )
        {
            idType = "1";
        }
        else  if(firstchar == "g" )
        {
            idType = "2";
        }
        else  if(firstchar == "e" )
        {
            idType = "3";
        }
    }
    if (isDebug == "true")
    {
        print  "idTypeList=",idTypeList, "idType=",idType, "id=",$1
    }
    if (match(idTypeList,idType) > 0)
    {
        print $1
    }
}
END{
    }
