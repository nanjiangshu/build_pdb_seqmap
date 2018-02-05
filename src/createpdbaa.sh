#!/bin/bash
# create the pdbaa files from SEQRES records in PDB files
# the default output path is $DATADIR/pdbaa
# the chain ID naming rule is the same as in file pdb_seqres.txt
# that is PDBID_chainID PDBID is in lowercase while chainID is in the original
# form as in the PDB file
# 2008-02-04, Nanjiang

# Copyright© Nanjiang Shu 
# Department of Materials and Environmental Chemistry, Stockholm University,
# Sweden
# Email: nanjiang.shu@mmk.su.se

function PrintHelp()
{
    echo "Usage: createpdbaa.sh [options] idListFile"
    echo " idListFile contains a list of standard chain IDs"
    echo "options:"
    echo "  -d path  : set the output path, default = $DATADIR/pdbaa"
    echo "  -n       : do not check the existence of necessary programs"
    echo "  -k|--keep: do not delete the temporary files"
    echo "  -h|--help: print this help message and exit"
    echo
    echo "Created 2007-08-11, last modified 2010-05-27 , Nanjiang"
}

if [ $# -lt 1 ]; then
    PrintHelp
    exit
fi

outpath=$DATADIR/pdbaa
isCheckExistence=true
isClean=true

while [ "$1" != "" ]; do
    case $1 in
        -h|--help)PrintHelp;exit;;
        -k|--keep)isClean=false;;
        -n)isCheckExistence=false;;
        -d)outpath=$2;shift;;
        *) idListFile=$1;;
    esac
    shift
done

if [ ! -f $idListFile ]; then
    echo "Error! Can not open file $idListFile"
    exit
fi

# check the necessary programs
if [ "$isCheckExistence" == "true" ]; then
    if [ "`which $BINPATH/lowercase > /dev/null`" != "" ]; then
        echo "Error! program \'lowercase\' does not exist"
        exit
    fi

    if [ "`which $BINPATH/splitfasta > /dev/null`" != "" ]; then
        echo "Error! program \'splitfasta\' does not exist"
        exit
    fi

    if [ "`which $BINPATH/getseqresseq > /dev/null`" != "" ]; then
        echo "Error! program \'getseqresseq\' does not exist"
        exit
    fi
fi

errFile=.$idListFile.$$.err.tmp
pdbfilelist=.$idListFile.$$.pdbfilelist.tmp
chainIDList=.$idListFile.$$.chainIDList.tmp
cPDBAAFile=.$idListFile.$$.cPDBAA.tmp
fastaFile=.$idListFile.$$.fasta.tmp

$BINPATH/getpdbfilepath -l $idListFile 1> $pdbfilelist 2> $errFile

if [  -s "$errFile" ]; then  #errFile exists and have a size greater than zero
    cat $errFile
    exit
fi

awk '{print substr($1,5,1)}' $idListFile > $chainIDList
paste $pdbfilelist $chainIDList > $cPDBAAFile
$BINPATH/getseqresseq -l $cPDBAAFile  > $fastaFile
$BINPATH/splitfasta -d $outpath $fastaFile

if [ "$isClean" == "true" ]; then 
    rm -f $errFile
    rm -f $pdbfilelist
    rm -f $chainIDList
    rm -f $cPDBAAFile
    rm -f $fastaFile
fi

#the following code is very slow
#((i=0))
#for id in $(cat $idListFile); do
#    chainID=${id:4:1}
#    pdbid=${id:0:4}
#    pdbid=`echo $pdbid | lowercase`
#    echo -e "$i \t getseqresseq `getpdbfilepath $id`  $chainID > $outpath/${pdbid}_${chainID}.aa"
#    getseqresseq `getpdbfilepath $id`  $chainID > $outpath/${pdbid}_${chainID}.aa
#    ((i++))
#done
