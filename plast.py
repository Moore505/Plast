#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import numpy as np

##################################################################################
#                                      Globals                                   #
##################################################################################
debug = False
MATCH_SCORE = 5
MISMATCH_PENALITY = -4
LAMBDA = 0.192
K = 0.176

##################################################################################
#                                      Classes                                   #
##################################################################################

# For a single sequence in DB, we have the genes and the corresponding details
class Sequence:

    # Constructor with name (>I|gat|...)
    def __init__(self, name, sequence):
        self.name = name.rstrip("\n")
        self.genes = sequence.rstrip("\n")
        self.kmers = getKmers(self.genes)

# We hold the database in here (database sequences)
class Database:

    # Constructor with db file name to initalise immediately
    def __init__(self, dbFile):
        self.sequences = []
        self.m = 0
        # Read and initiate sequences
        self.__init_database(dbFile)

    def __init_database(self, dbFile):
        if debug : print("(#1.1) Initialising database...")
        file = open(dbFile)
        lines = iter(file)
        if debug : print("(#1.2) Generating database k-mers...")
        for line in lines:
            currentSequence = Sequence(line, next(lines))
            self.sequences.append(currentSequence)
            # We keep a counter of every single gene for evalue
            self.m += len(currentSequence.genes)
        file.close

# We generate and compute every HSP in here
class HSPs:
    
    # Constructor with database, query and k (len of seed)
    def __init__(self, database, query, k):
        self.initial_HSPs = []
        self.extended_HSPs = []
        self.merged_HSPs = []
        self.final_HSPs = []
        self.database = database
        self.k = k # len(seed)
        self.query = query
        if debug : print("(#1.2) Generating query k-mers...")
        self.query_kmers = getKmers(self.query)
        self.__computeInitialHSPs(self.database.sequences, self.query_kmers)
        self.__extend_HSPs()
        self.__combine_extended_HSPs()
        self.__calculate_scores()

    # We first generate our initial HSP with perfect matches
    # between kmers from query and kmers from database sequences
    # (1.3 dans l'énoncé)
    def __computeInitialHSPs(self, db_sequences, query_kmers):
        if debug:
            print("(#1.3) Generating initial HSP...", end="\n\n")

        # for every query kmer, compare with each database kmer
        for q_kmer_index, q_kmer in enumerate(query_kmers):
            # We iterate trough database sequence
            for db_sequence_index, db_sequence in enumerate(db_sequences):
                # Then we iterate trough each kmer in current sequence
                # If we have an exact match, we save the informations
                # ** Custom seed not implemented yet **
                for db_sequence_kmer_index, db_kmer in enumerate(db_sequence.kmers):
                    if db_kmer == q_kmer:
                        # Since our seed is '11111111111', score is automatically 5 * k (for now)
                        initial_score = 5 * self.k
                        self.initial_HSPs.append([
                            q_kmer_index,            # Index of kmer in the list of query kmers (equal to starting index in original query)
                            db_sequence_index,       # Index of db_sequence in the list of db_sequences
                            db_sequence_kmer_index,  # Index of kmer in our list of kmers in db_sequence (equal to starting index in original sequence)
                            initial_score            # Initial score of the HSP
                            ])

                        if debug:
                            print("-" * 15 + " Match found " + "-" * 15)
                            print("q_kmer_index : " + str(q_kmer_index))
                            print("db_sequence_index : " + str(db_sequence_index))
                            print("db_sequence_kmer_index : " + str(db_sequence_kmer_index))
                            print(db_kmer)
                            print(q_kmer, end="\n\n")

    # We extend our initial HSPs here
    # (1.4 dans l'énoncé)
    def __extend_HSPs(self):
        if debug:
            print("(#1.4) Computing extended HSPs...", end="\n\n")

        # For every match, we extend the 2 kmers simultaneously to calculate the score
        for current_hsp_index, current_hsp in enumerate(self.initial_HSPs):

            # For readability
            q_kmer_index = current_hsp[0]
            db_sequence_index = current_hsp[1]
            db_sequence_kmer_index = current_hsp[2]
            initial_score = current_hsp[3]
            sequence_genes = self.database.sequences[db_sequence_index].genes # Complete db_sequence genes

            if debug:
                # Debugging current HSP informations
                print(("-" * 17) + " Match " + str(current_hsp_index+1) + " " + ("-" * 17))
                print("Initial score : " + str(initial_score))
                print("Original query and corresponding kmer : ")
                print(self.query)
                print(" " * q_kmer_index + self.query_kmers[q_kmer_index], end="\n\n")
                print("Original database sequence and corresponding kmer : ")
                print(self.database.sequences[db_sequence_index].name)
                print(sequence_genes)
                print(" " * db_sequence_kmer_index + self.database.sequences[db_sequence_index].kmers[db_sequence_kmer_index], end="\n\n")

            # We memorize the score
            highest_score = initial_score
            current_score = highest_score

            # Indexes to avoid OOB
            q_left = q_kmer_index - 1
            q_right = q_kmer_index + self.k + 1
            db_left = db_sequence_kmer_index - 1
            db_right = db_sequence_kmer_index + self.k + 1

            canGoLeft = (q_left > 0) and (db_left > 0)
            canGoRight = (q_right < len(self.query)) and (db_right < len(sequence_genes))

            while(canGoLeft or canGoRight):

                # Needs to be defined if kmer is already at most left of query/sequence
                go_left_score = 0

                # Left extension
                if canGoLeft:                    
                    # If next left character in query and sequence match
                    if self.query[q_left] == sequence_genes[db_left]:
                        go_left_score = current_score + MATCH_SCORE
                    else:
                        go_left_score = current_score + MISMATCH_PENALITY

                # Avoid infinite loop
                else: go_left_score = 0

                # Right extension
                if canGoRight:
                    # If next right character in query and sequence match
                    if self.query[q_right] == sequence_genes[db_right]:
                        go_right_score = current_score + MATCH_SCORE
                    else:
                        go_right_score = current_score + MISMATCH_PENALITY

                # Avoid infinite loop
                else: go_right_score = 0

                # We check if score downgraded beyond limit (highest - current > E)
                # We don't save the last extension in this case
                # Otherwise, we adjust the indexes and current score
                if go_left_score >= go_right_score:
                    if highest_score - go_left_score >= min_thres: break
                    q_left -= 1
                    db_left -= 1
                    current_score = go_left_score
                else:
                    if highest_score - go_right_score >= min_thres: break
                    q_right += 1
                    db_right += 1
                    current_score = go_right_score

                # We check if current score upgraded
                if current_score > highest_score:
                    highest_score = current_score

                # We re-evaluate if we can go right or left
                canGoLeft = (q_left > 0) and (db_left > 0)
                canGoRight = (q_right < len(self.query)) and (db_right < len(sequence_genes))

            # Once we stop, our indexes need to be adjusted
            # It will be a mismatch at the beginning and at the end
            q_left += 1
            db_left += 1
            q_right -=1
            db_right -=1

            if debug:
                print("Score : " + str(current_score))
                print(str(q_left) + " " + self.query[q_left:q_right] + " " + str(q_right))
                print(str(db_left) + " " + sequence_genes[db_left:db_right] + " " + str(db_right), end="\n\n")

            # Now we have indexes of extended HSP
            # Score will be recalculated after the merging
            self.extended_HSPs.append([
                q_left,
                q_right,
                db_left,
                db_right,          
                db_sequence_index,       # Index of db_sequence in the list of db_sequences
            ])          

    # We combine our extended HSPs here
    # (1.5 dans l'énoncé)
    def __combine_extended_HSPs(self):

        if debug: print("(#1.5) Combining extended HSPs...")

        # To easily split by db_sequence_index
        extended_hsps = np.array(self.extended_HSPs)
        extended_hsps = [extended_hsps[extended_hsps[:,4]==k] for k in np.unique(extended_hsps[:,4])]
        
        # Then we sort by start index for each db_sequence and merge overlapping indices
        for hsps_by_db_sequence in extended_hsps:
            sorted = hsps_by_db_sequence[hsps_by_db_sequence[:, 0].astype(int).argsort()]
            merged = [sorted[0]]
            for current in sorted:
                previous = merged[-1]
                if current[0] <= previous[1]:
                    previous[1] = max(previous[1], current[1])
                else:
                    merged.append(current)

            # merged will now contain the merged HSPs of the current sequence
            self.merged_HSPs.append(merged)

    # We calculate our scores on our final HSPs
    # (#1.6 dans l'énoncé)
    def __calculate_scores(self):
        if debug : print("(#1.6) Calculating scores for final HSPs...")
        final_unsorted = []
        for hsps_by_db_sequence in self.merged_HSPs:
            current_sequence = []
            for hsp in hsps_by_db_sequence:
                query = self.query[hsp[0]:hsp[1]]
                sequence = self.database.sequences[hsp[4]]
                sequence_aligned = sequence.genes[hsp[2]:hsp[2]+len(query)]
                score = HSP_score(query, sequence_aligned)
                bitscore = HSP_bitscore(score)
                evalue = HSP_evalue(self.database.m, len(self.query), bitscore)
                # We keep only if evalue < ss
                if evalue < signification_thres:
                    current_sequence.append([
                        evalue,
                        bitscore,
                        score,
                        hsp[0],
                        hsp[1],
                        query,
                        hsp[2],
                        hsp[3],
                        sequence.name,
                        sequence_aligned
                        ])
                    
            if len(current_sequence) > 0:
                x = np.array(current_sequence)
                sorted_array = x[x[:, 0].astype(np.float64).argsort()]
                final_unsorted.append(sorted_array[0])

        # Then we sort the best of each sequences
        if len(final_unsorted) > 0:
            x = np.array(final_unsorted)
            self.final_HSPs = x[x[:, 0].astype(np.float64).argsort()]

    def __str__(self) -> str:
        ret = ""
        # final_HSPs contains an array for each sequence HSPs
        for hsps_by_sequence in self.final_HSPs:
            # They are sorted by best evalue, so the first one is the best
            ret += hsps_by_sequence[8] + "\n"
            ret += "# Best HSP score:" + str(hsps_by_sequence[2]) + ", "
            ret += "bitscore:" + str(hsps_by_sequence[1]) + ", "
            ret += "evalue:" + str(hsps_by_sequence[0]) + "\n"

            ret += str(hsps_by_sequence[3]) + " "         # q_start index
            ret += str(hsps_by_sequence[5]) + " "         # query
            ret += str(hsps_by_sequence[4]) + "\n"        # q_end index
            ret += str(hsps_by_sequence[6]) + " "         # db_start index
            ret += str(hsps_by_sequence[9]) + " "         # sequence
            ret += str(hsps_by_sequence[7]) + "\n\n"      # db_end index
        ret += ("-" * 50) + "\n"
        ret += "Total : " + str(len(self.final_HSPs))
        return ret


##################################################################################
#                                Global functions                                #
##################################################################################

def parseArguments(debug=False):
    # Instantiate the parser
    parser = argparse.ArgumentParser(description='Plast : Primitive Local Alignment Search Tool')
    parser.add_argument('-db', type=str, required=True, help='Database used for search (required)')
    parser.add_argument('-i', type=str, required=True, help='Sequence query (required)')
    parser.add_argument('-E', type=int, help='Threshold (optionnal)')
    parser.add_argument('-ss', type=int, help='Signification threshold (optionnal)')
    parser.add_argument('-seed', type=str, help='Seed (optionnal)')
    args = parser.parse_args()

    # Declaring globals for accessing (to not have None if using defaults)
    global database_file
    database_file = ""
    global query
    query = ""
    global min_thres
    min_thres = 4
    global signification_thres
    signification_thres = 0.001
    global seed
    seed = '11111111111'

    # Assign variables and print if debug
    if debug : print('Database file: ' + args.db)
    database_file = args.db

    if debug : print('Sequence query: ' + args.i)
    query = args.i

    if args.E is not None:
        if debug : print('Treshold: ' + str(args.E))
        min_thres = args.E

    if args.ss is not None:
        if debug : print('Signification threshold: ' + str(args.ss))
        signification_thres = args.ss

    if args.seed is not None:
        if debug : print('Seed: ' + args.seed + '\n')
        seed = args.seed


# Generates kmers with k = len(seed)
# (1.2 dans l'énoncé)
def getKmers(sequence):
    kmers = []
    for i in range(0, len(sequence)):
	        if i+len(seed) < len(sequence):
		        kmers.append(sequence[i:i+len(seed)])
    return kmers

def HSP_score(s1, s2):
    score = 0
    if len(s1) == len(s2):
        for i in range(len(s1)):
            if s1[i] == s2[i] : score += MATCH_SCORE
            else : score += MISMATCH_PENALITY
    return score

def HSP_bitscore(score): return round((LAMBDA*score - np.log(K))/np.log(2))

# m : length of complete database genes
# n : length of complete query
def HSP_evalue(m, n, bitscore): return m*n*pow(2, -bitscore)

    
##################################################################################
#                                      Main                                      #
##################################################################################

if __name__ == "__main__":
    print()
    parseArguments(False)
    database = Database(database_file)
    hsps = HSPs(database, query, len(seed))
    if debug: print("(#1.7) Output of our final HSPs...", end="\n\n")
    print(hsps)