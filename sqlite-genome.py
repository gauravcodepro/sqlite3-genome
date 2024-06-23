#! usr/bin/env python3
# Author Gaurav 
# Universitat Potsdam
# Date 2024-6-20

import logging
import os
import sqlite3

import pandas as pd
import polars as pl


def sqlitewrapper(pathfasta, databasename):
  """
  Creates a sqlite database from a fasta file
  and build a wrapper to the cratedb for the genome
  render
  """
  if pathfasta and databasename is not None:
    readfile = '/'.join([os.path, pathfasta])
    connection = sqlite3.connect(databasename)
    readgenome = [i.strip() for i in open(pathfasta, "r").readlines()]
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
  with open(os.path.join(os.getcwd(),pathgff)) as gffread:
    with open(os.path.join(os.getcwd(),pathgff, ".mod.gff")) as gffwrite:
      store_gff = [line for line in gffread if not line.startswith("#")]
      # add a subfunction for the write of the headers with the modified information
      # for the loading into the database so that the crate, javascript database 
      # and all can sync
      gffwrite.write(storegff)
      gffread.close()
      gffwrite.close()
  readgff = pd.read_csv(os.path.join(os.getcwd(),pathgff, ".mod.gff"), sep = "\t")
  


if __name__ == __main__:
  sqlitewrapper(pathfasta, databasename)