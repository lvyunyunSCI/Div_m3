#!/usr/bin/env bash
ref_abb=$1 # reference genome abbreviation such as <HSAP> for representation of Homo sapiens (human)
ref_chrs_fas=$2 # reference genome fasta or fa file (only chomosomes were used)
qry_abb=$3 # query genome abbreviation as similar as $1
qry_chrs_fas=$4 # query genome fasta or fa file (only chomosomes were used)
thread=$5  # mutiple subprocess

#1 step # divide the chromsomes into every one
seqkit split -f -i --by-id-prefix '' --out-dir ${ref_abb}_split $ref_chrs_fas
seqkit split -f -i --by-id-prefix '' --out-dir ${query_abb}_split $qry_chrs_fas

#2 step # build reference Mash database 
find $PWD/${ref_abb}_split/*.* > ${ref_abb}.chrList
mash sketch -p $thread -k 31 -s 5000000000 -l ${ref_abb}.chrList -o $ref_abb # this step will build a mash database named as <${ref_abb}.msh>
mash info ${ref_abb}.msh > ${ref_abb}.msh.infor # get mash database infor for each chormosome

#3 step # with comparing of each chromosomes within query genome to mash database builded by reference chromosomes
find $PWD/${query_abb}_split/*.* > ${query_abb}.chrList
mash dist ${ref_abb}.msh -p $thread -s 5000000000 -l ${query_abb}.chrList > ${ref_abb}_${query_abb}_mashDistance # will get the Mash distance between query chromsomes and reference chromosomes

#4 step # with perl line commond to deal with the relationships between reference genome and query genome
cat  ${ref_abb}_${query_abb}_mashDistance | perl -lane '$hash->{$F[0]}->{$F[1]}=$F[2];END{for $ref(sort { $a cmp $b} keys %{$hash}){$n=0;for $query(sort {$hash->{$ref}->{$a} <=> $hash->{$ref}->{$b}} keys %{$hash->{$ref}}){$n++;$v=$hash->{$ref}->{$query};$newref=(split(/\//,$ref))[-1];$newquery=(split(/\//,$query))[-1];$newref=~s/.fa//g;$newquery=~s/.fa//g;print "$newref\t$newquery\t$v" if($n<=2)}}}'|sort -k1.6,1.7n -k3n >  ${ref_abb}_${query_abb}_mashDistance.filter # this will get a filer with first colum represent reference chromosomes, second colum represent the query chromosomes
cat ${ref_abb}_${query_abb}_mashDistance.filter |perl -lane 'BEGIN{print "Rchr\tQchr\tsubg\tMashD"};$count{$F[0]}++;print "$F[0]\t$F[1]\tSG$count{$F[0]}\t$F[2]"' > ${ref_abb}_${query_abb}_mashDistance.filter.Gadd # this will add a subgenome assignment and head line

#5 step # plot the result# Using of R (ggplot2 package) to plot the result
Rscript PlotMashDistance.R ${ref_abb}_${query_abb}_mashDistance.filter.Gadd ${ref_abb}_${query_abb}_mashDistance.filter.Gadd.pdf # final result *.Gadd.pdf
