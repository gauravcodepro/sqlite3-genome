#! usr/bin/env python3
# Author Gaurav
# Date 2024-6-20
# Creates a sqlite database from a fasta file, gene annotations,
# po_anntoations, geneannoationfiles, chromsome assembly files, protein annotation files.

import logging
import os
import sqlite3
import pandas as pd


def sqlitewrapper(pathfasta, databasename):
    if pathfasta and databasename is not None:
        readfile = '/'.join([os.path, pathfasta])
        connection = sqlite3.connect(databasename)
        readgenome = [i.strip() for i in open(readfile, "r").readlines()]
        fastagenome = {}
        for i in readgenome:
            if i.startswith(">"):
                path = i.strip()
                if i not in fastagenome:
                    fastagenome[i] = ""
                continue
            fastagenome[path] += i.strip()
        fastaseq = list(fastagenome.values())
        fastanames = [i.replace(">", "")for i in (list(fastagenome.keys()))]
        lenfasta = []
        for i in range(len(fastaseq)):
            lenfasta.append(len(fastaseq[i]))
    assert len(fastaseq) == len(fastanames) == len(lenfasta)
    connection = con.cursor()
    connection.execute("CREATE TABLE ids(fastaids)")
    connection.execute("CREATE TABLE seq(fastasequences)")
    connection.execute("CREATE TABLE length(fastalength)")
    connection.execute("""
       INSERT INTO ids VALUES
         tuple(([fastanames[i] for i in range(len(fastanames))]))
         """)
    connection.commit()
    connection.execute("""
       INSERT INTO seq VALUES
         tuple(([fastaseq[i] for i in range(len(fastaseq))]))
         """)
    connection.commit()
    connection.execute("""
       INSERT INTO length VALUES
         tuple(([lenfasta[i] for i in range(len(lenfasta))]))
         """)
    connection.commit()

    def dbprepare(gfffile, prepare):
        if gfffile and prepare == "yes":
            with open(gfffile, "r") as gffread:
                with open("gfffilemod", "w") as gffwrite:
                    gffwrite.write("column1" + "\t" + "column2" + "\t" +"column3" + "\t" +
                    "column4" + "\t" + "column5" + "\t" + "column6" + "\t" + "column7" +
                    "\t" + "column8" + "\t" + "column9\n")
                for line in gffread.readlines():
                    gffwrite.write(line)
                gffwrite.close()
        gffdataframe = pd.read_csv("gfffilemod", sep = "\t")
        genomeloc = gffdataframe["column1"].to_list()
        accession = gffdataframe["column2"].to_list()
        typeaccession = gffdataframe["column3"].to_list()
        start = gffdataframe["column4"].to_list()
        end = gffdataframe["column5"].to_list()
        length = gffdataframe["column6"].to_list()
        posstrand = gffdataframe["column7"].to_list()
        negstrand = gffdataframe["column8"].to_list()
        idlocation = gffdataframe["column9"].to_list()
        dbprepare = []
        for i in range(len(genomeloc)):
            dbprepare[genomeloc[i]] = [{genomeloc[i]: accession[i]}, {genomeloc[i]: typeaccession[i]},
                                        {genomeloc[i]: start[i]}, {genomeloc[i]:end[i]},
                                        {genomeloc[i]:length[i]}, {genomeloc[i]: posstrand[i]},
                                        {genomeloc[i]:negstrand[i]}, {genomeloc[i]: idlocation[i]}]dlocation[i]}])
        connection = con.cursor()
        connection.execute("CREATE TABLE genomedetails(genomeloc)")
        connection.execute("CREATE TABLE genomedetails(accession)")
        connection.execute("CREATE TABLE genomedetails(typeaccession)")
        connection.execute("CREATE TABLE genomedetails(start)")
        connection.execute("CREATE TABLE genomedetails(end)")
        connection.execute("CREATE TABLE genomedetails(length)")
        connection.execute("CREATE TABLE genomedetails(posstrand)")
        connection.execute("CREATE TABLE genomedetails(negstrand)")
        connection.execute("CREATE TABLE genomedetails(idlocation)")
        connection.execute("""
            INSERT INTO genomeloc VALUES
                     tuple(([genomeloc[i] for i in range(len(genomeloc))]))
         """)
        connection.commit()
         connection.execute("""
            INSERT INTO accession VALUES
                     tuple(([accession[i] for i in range(len(accession))]))
         """)
        connection.commit()
         connection.execute("""
            INSERT INTO typeaccession VALUES
                     tuple(([typeaccession[i] for i in range(len(typeaccession))]))
         """)
        connection.commit()
         connection.execute("""
            INSERT INTO start VALUES
                     tuple(([start[i] for i in range(len(start))]))
         """)
        connection.commit()
         connection.execute("""
            INSERT INTO end VALUES
                     tuple(([end[i] for i in range(len(end))]))
         """)
        connection.commit()
         connection.execute("""
            INSERT INTO length VALUES
                     tuple(([length[i] for i in range(len(length))]))
         """)
        connection.commit()
         connection.execute("""
            INSERT INTO poststrand VALUES
                     tuple(([posstrand[i] for i in range(len(posstrand))]))
         """)
        connection.commit()
         connection.execute("""
            INSERT INTO negsstrand VALUES
                     tuple(([negsstrand[i] for i in range(len(negsstrand))]))
         """)
        connection.commit()
        connection.execute("""
            INSERT INTO idlocation VALUES
                     tuple(([idlocation[i] for i in range(len(idlocation))]))
         """)
        connection.commit()


    def exonparse(pathgff):
        if pathgff:
            with open(gfffile, "r") as gffread:
                with open("gfffilemod", "w") as gffwrite:
                    gffwrite.write("column1" + "\t" + "column2" + "\t" +"column3" + "\t" +
                  "column4" + "\t" + "column5" + "\t" + "column6" + "\t" + "column7" +
                  "\t" + "column8" + "\t" + "column9\n")
                    for line in gffread.readlines():
                        gffwrite.write(line)
                    gffwrite.close()
            gffdataframe = pd.read_csv("gfffilemod", sep = "\t")
            select1 = gffdataframe["column3"].to_list()
            select2 = gffdataframe["column4"].to_list()
            select3 = gffdataframe["column5"].to_list()
            prepareexon = {}
            for i in range(len(select1)):
                if select1 == "exon":
                    prepareexon[select1[i]] = [{select1[i]: select2[i]}, {select1[i]: select2[i]}]
        connection = con.cursor()
        connection.execute("CREATE TABLE exon(select1)")
        connection.execute("CREATE TABLE exon(select2)")
        connection.execute("CREATE TABLE exon(select3)")
        connection.execute("""
            INSERT INTO select1 VALUES
                     tuple(([select1[i] for i in range(len(select1))]))
         """)
        connection.commit()
         connection.execute("""
            INSERT INTO select2 VALUES
                     tuple(([select2[i] for i in range(len(select2))]))
         """)
        connection.commit()
        connection.execute("""
            INSERT INTO select3 VALUES
                     tuple(([select3[i] for i in range(len(select3))]))
         """)
        connection.commit()

    def intronparse(pathgff):
        if pathgff:
            with open(gfffile, "r") as gffread:
                with open("gfffilemod", "w") as gffwrite:
                    gffwrite.write("column1" + "\t" + "column2" + "\t" +"column3" + "\t" +
                  "column4" + "\t" + "column5" + "\t" + "column6" + "\t" + "column7" +
                  "\t" + "column8" + "\t" + "column9\n")
                    for line in gffread.readlines():
                        gffwrite.write(line)
                    gffwrite.close()
            gffdataframe = pd.read_csv("gfffilemod", sep = "\t")
            select1 = gffdataframe["column3"].to_list()
            select2 = gffdataframe["column4"].to_list()
            select3 = gffdataframe["column5"].to_list()
            prepareintron = {}
            for i in range(len(select1)):
                if select1 == "exon":
                    prepareintron[select1[i]] = [{select1[i]: select2[i]}, {select1[i]: select3[i]}]
        connection = con.cursor()
        connection.execute("CREATE TABLE intron(select1)")
        connection.execute("CREATE TABLE intron(select2)")
        connection.execute("CREATE TABLE intron(select3)")
        connection.execute("""
            INSERT INTO select1 VALUES
                     tuple(([select1[i] for i in range(len(select1))]))
         """)
        connection.commit()
         connection.execute("""
            INSERT INTO select2 VALUES
                     tuple(([select2[i] for i in range(len(select2))]))
         """)
        connection.commit()
        connection.execute("""
            INSERT INTO select3 VALUES
                     tuple(([select3[i] for i in range(len(select3))]))
         """)
        connection.commit()

    def exonseq(pathgff, pathfasta):
        if pathgff and pathfasta:
            readfasta = [i.strip() for i in open(pathfasta, "r").readlines()]
            fastaseq = {}
            for i in readfasta:
                if i.startswith(">"):
                    path = i.strip()
                    if i not in fastaseq:
                        fastaseq[i] = ""
                    continue
                fastaseq[path] += i.strip()
            fasta_seq = list(fastaseq.values())
            fasta_names = [i.replace(">", "")for i in (list(fastaseq.keys()))]
            fastaparsedict = {}
            for i in range(len(fasta_seq)):
                fastaparsedict[fasta_seq[i]] = fasta_names[i]
            with open(gfffile, "r") as gffread:
                with open("gfffilemod", "w") as gffwrite:
                    gffwrite.write("column1" + "\t" + "column2" + "\t" +"column3" + "\t" +
                  "column4" + "\t" + "column5" + "\t" + "column6" + "\t" + "column7" +
                  "\t" + "column8" + "\t" + "column9\n")
                    for line in gffread.readlines():
                        gffwrite.write(line)
                    gffwrite.close()
            gffdataframe = pd.read_csv("gfffilemod", sep = "\t")
            exonpresent = gffdataframe["column1"].to_list()
            select1 = gffdataframe["column3"].to_list()
            select2 = gffdataframe["column4"].to_list()
            select3 = gffdataframe["column5"].to_list()
            exonseq ={}
            for i in range(len(select1)):
                if select1 == "exon":
                    exonseq[exonpresent[i]] = [select1[i], select2[i], select3[i]]
            exonseqprepare = {}
            exonseqkeys = list(exonseq.keys())
            exonseqvalues = list(exonseq.values())
            for i in range(len(exonseqkeys)):
                for j in range(len(fasta_seq)):
                    if exonseqkeys[i] == fasta_names[i]:
                        exonseqprepare[exonseqkeys[i]] == fasta_seq[i][exonseqvalues[i][0]:exonseqvalues[i][1]]
        connection = con.cursor()
        connection.execute("CREATE TABLE exonseq(fastaids)")
        connection.execute("CREATE TABLE exonseq(fastaseq)")
        connection.execute("""
            INSERT INTO fastaids VALUES
                     tuple(([exonseqprepare.keys()[i] for i in range(len(list(exonseqprepare.keys()))))]))
         """)
        connection.commit()
        connection.execute("""
            INSERT INTO select2 VALUES
                     tuple(([exonseqprepare.values()[i] for i in range(len(list(exonseqprepare.values()))))]))
         """)
        connection.commit()


    def intronseq(pathgff, pathfasta) :
        if pathgff and pathfasta:
            readfasta = [i.strip() for i in open(pathfasta, "r").readlines()]
            fastaseq = {}
            for i in readfasta:
                if i.startswith(">"):
                    path = i.strip()
                    if i not in fastaseq:
                        fastaseq[i] = ""
                    continue
                fastaseq[path] += i.strip()
            fasta_seq = list(fastaseq.values())
            fasta_names = [i.replace(">", "")for i in (list(fastaseq.keys()))]
            fastaparsedict = {}
            for i in range(len(fasta_seq)):
                fastaparsedict[fasta_seq[i]] = fasta_names[i]
            with open(gfffile, "r") as gffread:
                with open("gfffilemod", "w") as gffwrite:
                    gffwrite.write("column1" + "\t" + "column2" + "\t" +"column3" + "\t" +
                  "column4" + "\t" + "column5" + "\t" + "column6" + "\t" + "column7" +
                  "\t" + "column8" + "\t" + "column9\n")
                    for line in gffread.readlines():
                        gffwrite.write(line)
                    gffwrite.close()
            gffdataframe = pd.read_csv("gfffilemod", sep = "\t")
            intronpresent = gffdataframe["column1"].to_list()
            select1 = gffdataframe["column3"].to_list()
            select2 = gffdataframe["column4"].to_list()
            select3 = gffdataframe["column5"].to_list()
            intronseq ={}
            for i in range(len(select1)):
                if select1 == "exon":
                    intronseq[intronpresent[i]] = [select1[i], select2[i], select3[i]]
            intronseqprepare = {}
            intronseqkeys = list(intronseq.keys())
            intronseqvalues = list(intronseq.values())
            for i in range(len(intronseqkeys)):
                for j in range(len(fasta_seq)):
                    if intronseqkeys[i] == fasta_names[i]:
                        intronseqprepare[intronseqkeys[i]] == fasta_seq[i][intronseqvalues[i][0]:intronseqvalues[i][1]]
        connection = con.cursor()
        connection.execute("CREATE TABLE intronseq(fastaids)")
        connection.execute("CREATE TABLE intronseq(fastaseq)")
        connection.execute("""
            INSERT INTO fastaids VALUES
                     tuple(([intronseqprepare.keys()[i] for i in range(len(list(intronseqprepare.keys()))))]))
         """)
        connection.commit()
        connection.execute("""
            INSERT INTO select2 VALUES
                     tuple(([intronseqprepare.values()[i] for i in range(len(list(intronseqprepare.values()))))]))
         """)
        connection.commit()
                

    def goparsemongo(go_anntoations):
        with open(go_anntoations, "r") as goannotate:
            with open("goannotatemod.txt", "w") as gowrite:
                for line in goannotate.readlines():
                    if "!" not in line:
                        gowrite.write(line)
                    goannotatewrite.close()
        pareontologycomb1 = []
        prepareontologycomb2 = []
        prepareontologycomb3 = []
        repareontologycomb4 = []
        prepareontologycomb5 = []
        with open("goannotatemod.txt", "r") as gowrite:
            for line in gowrite.readlines():
                if len(line) <= 3:
                    continue
                else:
                    prepareontologycomb1.append({line.strip().split("\t")[0]: line.strip().split("\t")[1]})
                    prepareontologycomb2.append({line.strip().split("\t")[0]: line.strip().split("\t")[2]})
                    prepareontologycomb3.append({line.strip().split("\t")[0]: line.strip().split("\t")[3]})
                    prepareontologycomb4.append({line.strip().split("\t")[1]: line.strip().split("\t")[2]})
                    prepareontologycomb5.append({line.strip().split("\t")[1]: line.strip().split("\t")[3]})
