#!/bin/bash
# create dssp file given a list of pdbids or pdbfiles

# Copyright© Nanjiang Shu 
# Department of Materials and Environmental Chemistry, Stockholm University,
# Sweden
# Email: nanjiang.shu@mmk.su.se

function PrintHelp()
{
    echo "Usage: createdssp [options] pdbfile"
    echo " options"
    echo "  -d            : outpath, default = ./"
    echo "  --idlist file : pdbid list file"
    echo "  -l       file : pdbfile list file"
    echo " --help|-h      : print this help message and exit"
    echo
    echo " Created 2007-10-21, updated 2010-05-27, Nanjiang Shu"
}

function PDBFileName2PDBID()# $pdbfile#{{{
{
    local pdbfile=$1
    basename=`basename "$pdbfile"`
    rootname=${basename%.*}
    pdbid=${rootname#pdb} 
    echo $pdbid | tr "[:upper:]" "[:lower:]"
}
#}}}
function CreateDSSP() # pdbfile num pdbid#{{{
{
    local pdbfile=$1
    local num=$2
    local pdbid=`echo $3 | tr "[:upper:]" "[:lower:]"` 
    if [ "$pdbid" == "" ]; then
        pdbid=`PDBFileName2PDBID $pdbfile`
    fi 
    outfile=$outpath/$pdbid.dssp
    $BINPATH/dssp -ssa $pdbfile $outpath/$pdbid.dssp
    echo -e "$num \t the DSSP file output to $outfile/"
}
#}}}
if [ $# -lt 1 ]; then
    PrintHelp
    exit
fi

outpath=./
pdbidListFile=
pdbFileListFile=
pdbFileList=

while [ "$1" != "" ]; do
    case $1 in
        --help|-h) PrintHelp; exit;;
        -d) outpath=$2; shift;;
        --idlist) pdbidListFile=$2; shift;;
        -l) pdbFileListFile=$2; shift;;
        *) pdbFileList="$pdbFileList $1";;
    esac
    shift
done

# check the existance of the program dssp
prog=$BINPATH/dssp
if ! which $prog > /dev/null; then
    echo "Error! the program $prog does not exist, exit"  >&2
    exit
fi
prog=$BINPATH/getpdbfilepath
if ! which $prog > /dev/null; then
    echo "Error! the program $prog does not exist, exit"  >&2
    exit
fi

if [ ! -d $outpath ]; then
    mkdir -p $outpath
fi

if [ "$pdbFileList" != "" ]; then 
    ((i=0))
    for pdbfile in $pdbFileList; do
        CreateDSSP $pdbfile $i 
        ((i++))
    done
elif [ "$pdbidListFile" != "" ]; then
    if [ ! -f "$pdbidListFile" ]; then 
        echo "Error! The pdbidListFile $pdbidListFile does not exit." >&2
        exit
    fi
    ((i=0))
    for pdbid in $(cat $pdbidListFile); do
        pdbfile=`$BINPATH/getpdbfilepath $pdbid`
        if [[ "$pdbfile" =~ "no" ]]; then
            continue
        fi
        CreateDSSP $pdbfile $i $pdbid
        ((i++))
    done
elif [ "$pdbFileListFile" != "" ]; then
    if [ ! -f "$pdbFileListFile" ]; then 
        echo "Error! The pdbFileListFile $pdbFileListFile does not exit." >&2
        exit
    fi
    ((i=0))
    for pdbfile in $(cat $pdbFileListFile); do
        CreateDSSP $pdbfile $i
        ((i++))
    done
else
    echo "Error! supply pdbfile or pdbfilelist or pdbidlist"
fi
