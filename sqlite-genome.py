#! usr/bin/env python3
# Author Gaurav
# Universitat Potsdam
# Date 2024-6-20
# Creates a sqlite database from a fasta file, gene annotations, po_anntoations, geneannoationfiles, chromsome assembly files,
# protein annotation files.

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


def gffload(pathfasta, pathgff):
    with open(os.path.join(os.getcwd(), pathgff)) as gffread:
        with open(os.path.join(os.getcwd(), pathgff, ".mod.gff")) as gffwrite:
            store_gff = [line for line in gffread if not line.startswith("#")]
            gffwrite.write(store_gff)
            gffread.close()
            gffwrite.close()
    readgff = pd.read_csv(os.path.join(
        os.getcwd(), pathgff, ".mod.gff"), sep="\t")
