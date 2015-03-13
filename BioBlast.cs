/*
"The contents of this file are subject to the Mozilla Public License
Version 1.1 (the "License"); you may not use this file except in
compliance with the License. You may obtain a copy of the License at
http://www.mozilla.org/MPL/

Software distributed under the License is distributed on an "AS IS"
basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
License for the specific language governing rights and limitations
under the License.

The Original Code is the content of this file.

The Initial Developer of the Original Code is Dr. Robert Kofler.
Carrer del Baluard; 73 2°1a
08003 Barcelona
All Rights Reserved.

Contributor(s): Dr. Robert Kofler.

Alternatively, the contents of this file may be used under the terms
of the _____ license (the  "[___] License"), in which case the
provisions of [______] License are applicable instead of those
above. If you wish to allow use of your version of this file only
under the terms of the [____] License and not to allow others to use
your version of this file under the MPL, indicate your decision by
deleting the provisions above and replace them with the notice and
other provisions required by the [___] License. If you do not delete
the provisions above, a recipient may use your version of this file
under either the MPL or the [___] License."
*/
using System;
using System.Collections.Generic;
using System.Text;


namespace Bio.Blast
{
    using Bio.Seq;
    using Bio.Alignment;
    using Bio.Seq.IO;
    using Bio.Alignment.Misc;
    using Bio.Blast.Misc;

    [Flags]
    public enum BlastNOptions
    {
        None=0x00,
        PerformUngappedExtension=0x01,
        ConsiderSeedsOverlappingWithLastAlignment=0x02
    }
        /// <summary>
        /// Implementation of BLASTN
        /// Allows an extremely flexible design and custom tailoring of BLASTN.
        /// </summary>
    public class BlastN
    {
        //STRATEGIES
        private IBlastSeedProcessor seedProcessor;
        private ISignificator significator;
        private IBlastPostProcessor postproces;
        private IBlastHashtableConstructor hasher;

        //Progress Report EventHandler
        public event BioProgressReporter Reporter;
        private int progressToReport = 100;
        private int progressCounter = 0;

        //Prototype - Strategie for dynamic programming
        IAnchoredDynamicProgramming anchoredDynamic;

        //Options
        private bool ignoreSeedsOverlappingWithLastAlignment = true;
        private bool ungappedExtensionPhase = false;

        //HASHING
        private int wordSize;
        private DNADictionary<List<BlastHashSeed>> hashSeeds = null;

        //Working variables
        private SubstitutionMatrix substitutionMatrix;
        private List<NucleotideSequence> databaseSequences;


        private float minScoreForGapedExtension;
#if DEBUG
        private int countUngappedExtension = 0;
#endif




        public BlastN(List<NucleotideSequence> databaseSequences, SubstitutionMatrix substitutionMatrix)
        {
            this.databaseSequences = databaseSequences;
            this.substitutionMatrix = substitutionMatrix;
        }

        public BlastN(List<NucleotideSequence> databaseSequences, SubstitutionMatrix substitutionMatrix, BlastNOptions options)
        {
            this.databaseSequences = databaseSequences;
            this.substitutionMatrix = substitutionMatrix;

            if ((options & BlastNOptions.ConsiderSeedsOverlappingWithLastAlignment) == BlastNOptions.ConsiderSeedsOverlappingWithLastAlignment) this.ignoreSeedsOverlappingWithLastAlignment = false;
            if ((options & BlastNOptions.PerformUngappedExtension) == BlastNOptions.PerformUngappedExtension) this.ungappedExtensionPhase = true;
        }

        public BlastN(NucleotideSequence databaseSequence, SubstitutionMatrix substitutionMatrix)
        {
            this.databaseSequences = new List<NucleotideSequence>();
            this.databaseSequences.Add(databaseSequence);
            this.substitutionMatrix = substitutionMatrix;
        }

        /// <summary>
        /// Retrieve a List of all possible alignments with the given query sequence, 
        /// provided the hits fullfill the minimum requirements as specified in the ISignificator strategy
        /// </summary>
        /// <param name="query"></param>
        /// <returns></returns>
        public List<IPairwiseAlignmentContainer> GetAlignment(NucleotideSequence query)
        {
            //First initialize the dictionary
            if (hashSeeds == null) CreateHashTable();

            //Initialize the Dictionary and the word boundaries
            List<PairwiseNucleotideSequenceAlignment> toReturn = new List<PairwiseNucleotideSequenceAlignment>();


            //Plus Strand
            toReturn.AddRange(GetNucleotideSequenceAlignmentsSingleStrand(query, true));
            //Reverse Complement (Minus strand)
            toReturn.AddRange(GetNucleotideSequenceAlignmentsSingleStrand(query, false));
            
            //
            //Report progress
            //
            this.progressCounter++;
            if (Reporter != null && this.progressCounter % this.progressToReport == 0) Reporter(this, String.Format("Processed {0} sequences",progressCounter));
            
            //Postprocessor
            //Give the list of the high scoring pairs to the postprocessor and return the results to the user
            return this.postproces.GetAlignments(toReturn);
        }

        private void CreateHashTable()
        {
            if (Reporter != null) Reporter(this, String.Format("Creating a Hash table using a wordsize of {0} bp and a stepsize of {1} bp", hasher.Wordsize, hasher.Stepsize));

            if (this.databaseSequences.Count == 1) this.hashSeeds = hasher.GetBlastHashtable(databaseSequences[0]);
            else this.hashSeeds = hasher.GetBlastHashtable(databaseSequences);
        }


        private List<PairwiseNucleotideSequenceAlignment> GetNucleotideSequenceAlignmentsSingleStrand(INucleotideSequence query, bool plusStrand)
        {
            //The sequence which should be searched, the normal strand or the reverse complement
            ISequenceContainer querySequence;
            List<PairwiseNucleotideSequenceAlignment> alignments = new List<PairwiseNucleotideSequenceAlignment>();

            //Use the normal or the reverse complement strand
            if (plusStrand == true) querySequence = query.Sequence;
            else querySequence = NucleotideSequenceUtility.GetReverseComplement(query.Sequence);


            //Get the seeds
            List<BlastSeed> seeds = GetSeeds(querySequence);
            if (seeds == null) return null;
            seeds.Sort(new SortBlastSeed_DatabaseID_DatabasePosition_QueryPosition());

            int lastHitPosition = -1;
            int lastHitID = -1;


            //Test each seed
            for (int i = 0; i < seeds.Count; i++)
            {
                //
                //DECIDE whether or not the given seed should be tested either ungaped or dynamic
                //
                if (!ignoreSeedsOverlappingWithLastAlignment || seeds[i].DatabaseID != lastHitID || seeds[i].DatabaseStartPosition > lastHitPosition)
                {

                    NucleotideSequence databaseSequence = databaseSequences[seeds[i].DatabaseID];

                    //
                    //UNGAPED EXTENSION if desired
                    //
                    AnchoredUngappedExtension tester = null;
                    if (ungappedExtensionPhase == true)//Should an ungapped extension be performed
                    {
                        tester = new AnchoredUngappedExtension(databaseSequence.Sequence, querySequence, this.substitutionMatrix, seeds[i].DatabaseStartPosition + 1, seeds[i].QueryStartPosition + 1, 30);
#if DEBUG
                        countUngappedExtension += 1;
#endif
                    }

                    //
                    //Dynamic programming
                    //
                    if (tester == null || tester.Score >= minScoreForGapedExtension)
                    {

#if DEBUG
                        countUngappedExtension += 1;
#endif

                        IAnchoredDynamicProgramming smithWater = this.anchoredDynamic.NewDynamicProgramming(databaseSequence.Sequence, querySequence, this.substitutionMatrix, seeds[i].DatabaseStartPosition + 1, seeds[i].QueryStartPosition + 1, plusStrand);

                        PairwiseNucleotideSequenceAlignment al = new PairwiseNucleotideSequenceAlignment(smithWater.GetAlignment(), databaseSequence.Name, query.Name, smithWater.Start_Database, smithWater.Start_Query, smithWater.End_Database, smithWater.End_Query);
                        al.SubstitutionMatrix = this.substitutionMatrix;
                        al.Score = smithWater.Score;
                        al.PlusPlusStrand = plusStrand;
                        al.AlgorithmUsedForAlignment = smithWater.GetType();
                        al.Significator = this.significator;
                        al.LengthQueryParent = querySequence.Length;
                        al.LengthDatabaseParent = (int)databaseSequence.Length;

                        //If the alignment is significant using the given significator store it in the sequences which should be reported
                        if (al.IsSignificant())
                        {
                            alignments.Add(al);
                            lastHitID = seeds[i].DatabaseID;
                            lastHitPosition = al.EndDatabase - 1;
                        }
                    }




                }



            }


            return alignments;

        }




        /// <summary>
        /// Return a number of non overlapping seeds for the database sequence
        /// </summary>
        private List<BlastSeed> GetSeeds(ISequenceContainer query)
        {
            List<BlastSeed> protoSeeds = new List<BlastSeed>();


            //int lastSeedPosition = -1;

            //Wordnumber = n-w+1

            //01234567890123456789
            //AAAAAAAAAACCCCCCCCCC
            //the last has to be performed for <20-10+1 -> 10
            int range = query.Length - wordSize + 1;
            for (int i = 0; i < range; i++)
            {
                List<BlastHashSeed> prelim = null;
                ISequenceContainer sequencetotest = query.SubSequence(i, wordSize);
                if (hashSeeds.ContainsKey(sequencetotest))
                {
                    prelim = hashSeeds[sequencetotest];
                }

                if (prelim != null)
                {
                    //Write seeds into protoSeedList
                    foreach (BlastHashSeed k in prelim)
                    {
                        //k=databasePos; i=queryPos, both values are zero based array values
                        protoSeeds.Add(new BlastSeed(k, i));
                    }
                }
            }



            return this.seedProcessor.GetSeeds(protoSeeds);


        }



        /// <summary>
        /// Set or get the BlastSeedProcessor
        /// </summary>
        public IBlastSeedProcessor SeedProcessor
        {
            get
            {
                return this.seedProcessor;
            }
            set
            {
                this.seedProcessor = value;
            }
        }
        /// <summary>
        /// Get or set the hash-table constructing strategy
        /// </summary>
        public IBlastHashtableConstructor HashConstructor
        {
            get
            {
                return this.hasher;
            }
            set
            {
                this.hasher = value;
                this.wordSize = hasher.Wordsize;
            }
        }

        /// <summary>
        /// Get or set the postprocessor
        /// </summary>
        public IBlastPostProcessor PostProcessor
        {
            get
            {
                return this.postproces;
            }
            set
            {
                this.postproces = value;
            }
        }

        public ISignificator Significator
        {
            get
            {
                return this.significator;
            }
            set
            {
                this.significator = value;
            }
        }

        public IAnchoredDynamicProgramming AnchoredDynamicProgrammingAlgorithm
        {
            get
            {
                return this.anchoredDynamic;
            }
            set
            {
                this.anchoredDynamic = value;
            }

        }




        /// <summary>
        /// Minimum score for a ungaped extension to be passed to the dynamic programing algorithm
        /// </summary>
        public float MinimumScoreDuringGapedExtension
        {
            get
            {
                return this.minScoreForGapedExtension;
            }
            set
            {
                this.minScoreForGapedExtension = value;
            }
        }

        /// <summary>
        /// Get or set the progress which should be reported by the progress reporter
        /// </summary>
        public int ProgressToReport
        {
            get
            {
                return this.progressToReport;
            }
            set
            {
                this.progressToReport = value;
            }
        }


    }


        public class BlastNFactory
        {
            /// <summary>
            /// Obtain the default Blast implementation which can best be used for Gene-expression profiling using 454 sequencing
            /// </summary>
            /// <param name="databaseSequences">the list of the database sequences, in this case the genes</param>
            /// <returns></returns>
            public static BlastN GetExpressionProfiling454SequenceBlastN(List<NucleotideSequence> databaseSequences)
            {
                BlastN blast = new BlastN(databaseSequences, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(3.0F, 5.0F, 11.0F, 1.5F));
                blast.AnchoredDynamicProgrammingAlgorithm = new Anchored454SmithWatermanGotoh(40.0F);

                blast.Significator = new BlastSignificator_All();
                blast.HashConstructor = new BlastHashProcessor_NonOverlapping(11, 100);
                blast.SeedProcessor = new BlastSeedProcessor_BestDiagonal(3, 1, 5);
                blast.PostProcessor = new BlastPostProcessorHighestScore(1);
                return blast;

            }



            /// <summary>
            /// Obtain the default Blast implementation which can best be used for Gene-expression profiling using 454 sequencing
            /// </summary>
            /// <param name="databaseSequences">the list of the database sequences, in this case the genes</param>
            /// <returns></returns>
            public static BlastN GetExpressionProfiling454SequenceBlastN_Intron_KeepBestOverlap(List<NucleotideSequence> databaseSequences)
            {
                SubstitutionMatrix substitutionMatrix = SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(3.0F, 5.0F, 11.0F, 1.5F);
                BlastN blast = new BlastN(databaseSequences, substitutionMatrix);
                blast.AnchoredDynamicProgrammingAlgorithm = new Anchored454SmithWatermanGotoh(40.0F);
                blast.Significator = new BlastSignificator_All();
                blast.HashConstructor = new BlastHashProcessor_NonOverlapping(11, 100);
                blast.SeedProcessor = new BlastSeedProcessor_BestCumulativeDiagonal(2, 1, blast.HashConstructor, null, 5, 1);
                blast.PostProcessor = new BlastPostProcessor_PartialAlignmentAggregator(1, substitutionMatrix);
                return blast;

            }


            #region PanGEA related BlastN versions

            public static BlastN GetExpressionProfiling454SequenceBlastN_Intron_KeepBestOverlap(List<NucleotideSequence> databaseSequences, int wordsize, int minimumDiagonalLength, int lowComplexityThreshold, float hitScore, float mismatchPenalty, float gapExistPenalty, float gapExtendPenalty, float homopolymereTransgressionPenalty, bool useMaximumDistance, int maximumDistanceBetweenExons, float unambiguityScoreThreshold)
            {
                SubstitutionMatrix substitutionMatrix = SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(hitScore, mismatchPenalty, gapExistPenalty, gapExtendPenalty);
                BlastN blast = new BlastN(databaseSequences, substitutionMatrix);
                blast.AnchoredDynamicProgrammingAlgorithm = new Anchored454SmithWatermanGotoh(homopolymereTransgressionPenalty);

                blast.Significator = new BlastSignificator_All();
                blast.HashConstructor = new BlastHashProcessor_NonOverlapping(wordsize, lowComplexityThreshold);
                blast.ProgressToReport = 100;
                ///Should a maximum distance between exons be used or not
                if (useMaximumDistance)
                {
                    blast.SeedProcessor = new BlastSeedProcessor_BestCumulativeDiagonal(minimumDiagonalLength, 1, blast.HashConstructor, maximumDistanceBetweenExons, 5, 1);
                    blast.PostProcessor = new BlastPostProcessor_PartialAlignmentAggregator(1, substitutionMatrix, maximumDistanceBetweenExons, unambiguityScoreThreshold);
                }
                else //No maximum distance
                {
                    blast.SeedProcessor = new BlastSeedProcessor_BestCumulativeDiagonal(minimumDiagonalLength, 1, blast.HashConstructor, null, 5, 1);
                    blast.PostProcessor = new BlastPostProcessor_PartialAlignmentAggregator(1, substitutionMatrix, null, unambiguityScoreThreshold);
                }
                return blast;
            }


            public static BlastN GetExpressionProfiling454SequenceBlastN(List<NucleotideSequence> databaseSequences, int wordsize, int minimumDiagonalLength, int lowComplexityThreshold, float hitScore, float mismatchPenalty, float gapExistPenalty, float gapExtendPenalty, float homopolymereTransgressionPenalty, float unambiguityScoreThreshold)
            {
                BlastN blast = new BlastN(databaseSequences, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(hitScore, mismatchPenalty, gapExistPenalty, gapExtendPenalty));
                blast.AnchoredDynamicProgrammingAlgorithm = new Anchored454SmithWatermanGotoh(homopolymereTransgressionPenalty);
                blast.Significator = new BlastSignificator_All();
                blast.HashConstructor = new BlastHashProcessor_NonOverlapping(wordsize, lowComplexityThreshold);
                blast.SeedProcessor = new BlastSeedProcessor_BestDiagonal(minimumDiagonalLength, 1, 10);
                blast.PostProcessor = new BlastPostProcessorHighestScore(1, unambiguityScoreThreshold);
                blast.ProgressToReport = 100;
                return blast;
            }


            public static BlastN GetExpressionProfilingSequenceBlastN_Intron_KeepBestOverlap(List<NucleotideSequence> databaseSequences, int wordsize, int minimumDiagonalLength, int lowComplexityThreshold, float hitScore, float mismatchPenalty, float gapExistPenalty, float gapExtendPenalty, bool useMaximumDistance, int maximumDistanceBetweenExons, float unambiguityScoreThreshold)
            {
                SubstitutionMatrix substitutionMatrix = SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(hitScore, mismatchPenalty, gapExistPenalty, gapExtendPenalty);
                BlastN blast = new BlastN(databaseSequences, substitutionMatrix);
#if Flexible
                  blast.AnchoredDynamicProgrammingAlgorithm = new AnchoredDynamicBandedSmithWatermanGotoh(60,10,100);
#else
                blast.AnchoredDynamicProgrammingAlgorithm = new AnchoredSmithWatermanGotoh();
#endif
                blast.Significator = new BlastSignificator_All();
                blast.HashConstructor = new BlastHashProcessor_NonOverlapping(wordsize, lowComplexityThreshold);
                blast.ProgressToReport = 100;

                if (useMaximumDistance)
                {
                    blast.SeedProcessor = new BlastSeedProcessor_BestCumulativeDiagonal(minimumDiagonalLength, 1, blast.HashConstructor, maximumDistanceBetweenExons, 5, 1);
                    blast.PostProcessor = new BlastPostProcessor_PartialAlignmentAggregator(1, substitutionMatrix, maximumDistanceBetweenExons, unambiguityScoreThreshold);
                }
                else
                {
                    blast.SeedProcessor = new BlastSeedProcessor_BestCumulativeDiagonal(minimumDiagonalLength, 1, blast.HashConstructor, null, 5, 1);
                    blast.PostProcessor = new BlastPostProcessor_PartialAlignmentAggregator(1, substitutionMatrix, null, unambiguityScoreThreshold);
                }

                return blast;
            }


            /// <summary>
            /// Return an BlastN implementation
            /// </summary>
            /// <param name="databaseSequences"></param>
            /// <param name="wordsize"></param>
            /// <param name="minimumDiagonalLength"></param>
            /// <param name="lowComplexityThreshold"></param>
            /// <param name="hitScore"></param>
            /// <param name="mismatchPenalty"></param>
            /// <param name="gapExistPenalty"></param>
            /// <param name="gapExtendPenalty"></param>
            /// <returns></returns>
            public static BlastN GetExpressionProfilingSequenceBlastN(List<NucleotideSequence> databaseSequences, int wordsize, int minimumDiagonalLength, int lowComplexityThreshold, float hitScore, float mismatchPenalty, float gapExistPenalty, float gapExtendPenalty, float unambiguityScoreThreshold)
            {
                BlastN blast = new BlastN(databaseSequences, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(hitScore, mismatchPenalty, gapExistPenalty, gapExtendPenalty));
#if Flexible
                  blast.AnchoredDynamicProgrammingAlgorithm = new AnchoredDynamicBandedSmithWatermanGotoh(60,10,100);
#else
                blast.AnchoredDynamicProgrammingAlgorithm = new AnchoredSmithWatermanGotoh();
#endif
                blast.AnchoredDynamicProgrammingAlgorithm = new AnchoredSmithWatermanGotoh();
                blast.Significator = new BlastSignificator_All();
                blast.HashConstructor = new BlastHashProcessor_NonOverlapping(wordsize, lowComplexityThreshold);
                blast.SeedProcessor = new BlastSeedProcessor_BestDiagonal(minimumDiagonalLength, 1, 10);
                blast.PostProcessor = new BlastPostProcessorHighestScore(1, unambiguityScoreThreshold);
                blast.ProgressToReport = 100;
                return blast;
            }

            #endregion
        }




        /// <summary>
        /// Get a Blast hashtable, for each word the sequence id and the position within the sequence is stored.
        /// This implementation uses non overlapping words
        /// </summary>
        public class BlastHashProcessor_NonOverlapping : IBlastHashtableConstructor
        {
            private int wordsize;
            private int lowComplexityCutoff;
            /// <summary>
            /// Create a new Blast hash table constructor; non overlapping words will be used
            /// </summary>
            /// <param name="wordsize">the size of the word</param>
            public BlastHashProcessor_NonOverlapping(int wordsize, int lowComplexityCutoff)
            {
                this.wordsize = wordsize;
                this.lowComplexityCutoff = lowComplexityCutoff;
            }



            public DNADictionary<List<BlastHashSeed>> GetBlastHashtable(List<NucleotideSequence> databaseSequences)
            {
                DNADictionary<List<BlastHashSeed>> ret = new DNADictionary<List<BlastHashSeed>>(wordsize); //Create DNADictionary
                //
                ISequenceContainer template = SequenceFactory.GetDefaultSequence();
                if (databaseSequences.Count >= UInt16.MaxValue) throw new ArgumentOutOfRangeException(String.Format("Error in BlastHashProcessor.GetBlastHahstable(); the number of sequences exceeds the maximum allowed value {0}", Int16.MaxValue));

                long nucCount = 0;
                for (int i = 0; i < databaseSequences.Count; i++)
                {
                    ISequenceContainer data = databaseSequences[i].Sequence;
                    nucCount += data.Length;
                    for (int k = 0; k < data.Length - wordsize + 1; k += wordsize)
                    {
                        ISequenceContainer seq = data.ConvertTo(template,k, wordsize);
                        if (ret.ContainsKey(seq))
                        {
                            ret[seq].Add(new BlastHashSeed((ushort)i, k));
                        }
                        else
                        {
                            List<BlastHashSeed> list = new List<BlastHashSeed>();
                            list.Add(new BlastHashSeed((ushort)i, k));
                            ret.Add(seq, list);
                        }

                    }

                }


                //
                //Low complexity cutoff
                //
                int entries = (int)Math.Pow(4.0, this.wordsize);
                double expected = ((double)(nucCount / wordsize)) / (double)entries;
                //Expected is the expected number of occurences of a word with length wordsize


                int thresholdCount = (int)(expected * this.lowComplexityCutoff);
                if (thresholdCount < 5) thresholdCount = 5;

                for (int i = 0; i < entries; i++)
                {
                    List<BlastHashSeed> l = ret.GetValueForHash(i);
                    if (l != null)
                    {
                        if (ret.GetValueForHash(i).Count >= thresholdCount) ret[i] = null;
                        else l.TrimExcess(); //If the entry should be retained, trim excessively used capacity of the array!
                    }

                }


                return ret;


            }

            public DNADictionary<List<BlastHashSeed>> GetBlastHashtable(NucleotideSequence databaseSequence)
            {
                DNADictionary<List<BlastHashSeed>> ret = new DNADictionary<List<BlastHashSeed>>(wordsize);
                ISequenceContainer template = SequenceFactory.GetDefaultSequence();

                ISequenceContainer data = databaseSequence.Sequence;
                for (int k = 0; k < data.Length - wordsize + 1; k += wordsize)
                {
                    ISequenceContainer seq = data.ConvertTo(template,k, wordsize);
                    if (ret.ContainsKey(seq))
                    {
                        ret[seq].Add(new BlastHashSeed(k));
                    }
                    else
                    {
                        List<BlastHashSeed> list = new List<BlastHashSeed>();
                        list.Add(new BlastHashSeed(k));
                        ret.Add(seq, list);
                    }

                }

                //Remove low complexity
                int entries = (int)Math.Pow(4.0, this.wordsize);
                double expected = ((double)(data.Length / wordsize)) / (double)entries;
                //Expected is the expected number of occurences of a word with length wordsize


                int thresholdCount = (int)(expected * this.lowComplexityCutoff);
                if (thresholdCount < 5) thresholdCount = 5;

                for (int i = 0; i < entries; i++)
                {
                    List<BlastHashSeed> l = ret.GetValueForHash(i);
                    if (l != null)
                    {
                        if (ret.GetValueForHash(i).Count >= thresholdCount) ret[i] = null;
                        else l.TrimExcess(); //If the entry should be retained, trim excessively used capacity of the array!
                    }

                }



                return ret;
            }

            /// <summary>
            /// Retrieve the wordsize
            /// </summary>
            public int Wordsize
            {
                get
                {
                    return this.wordsize;
                }
            }
            public int Stepsize
            {
                get
                {
                    return this.wordsize;
                }
            }

            private DNADictionary<List<BlastHashSeed>> RemoveLowComplexity(DNADictionary<List<BlastHashSeed>> toRemove, long nucCount)
            {
                int entries = (int)Math.Pow(4.0, this.wordsize);
                double expected = ((double)(nucCount / wordsize)) / (double)entries;
                //Expected is the expected number of occurences of a word with length wordsize


                int thresholdCount = (int)(expected * this.lowComplexityCutoff);
                if (thresholdCount < 5) thresholdCount = 5;

                for (int i = 0; i < entries; i++)
                {
                    if (toRemove.GetValueForHash(i) != null && toRemove.GetValueForHash(i).Count >= thresholdCount)
                    {
                        toRemove[i] = null;
                    }
                }

                return toRemove;

            }
        }

        /// <summary>
        /// Get a Blast hashtable, for each word the sequence id and the position within the sequence is stored.
        /// This implementation uses overlapping words
        /// </summary>
        public class BlastHashProcessor_Overlapping : IBlastHashtableConstructor
        {

            private int wordsize;
            private int lowComplexityCutoff;
            public BlastHashProcessor_Overlapping(int wordsize, int lowComplexityCutoff)
            {
                this.wordsize = wordsize;
                this.lowComplexityCutoff = lowComplexityCutoff;
            }

            public DNADictionary<List<BlastHashSeed>> GetBlastHashtable(List<NucleotideSequence> databaseSequences)
            {
                DNADictionary<List<BlastHashSeed>> ret = new DNADictionary<List<BlastHashSeed>>(wordsize);
                ISequenceContainer template = SequenceFactory.GetDefaultSequence();
                long nucCount = 0;
                if (databaseSequences.Count >= ushort.MaxValue) throw new ArgumentOutOfRangeException(String.Format("Error in BlastUtility.GetDictionary(); the number of sequences exceeds the maximum allowed value {0}", ushort.MaxValue));


                for (int i = 0; i < databaseSequences.Count; i++)
                {
                    ISequenceContainer data = databaseSequences[i].Sequence;
                    nucCount += data.Length;
                    for (int k = 0; k < data.Length - wordsize + 1; k++)
                    {
                        ISequenceContainer seq = data.ConvertTo(template,k, wordsize);

                        if (ret.ContainsKey(seq))
                        {
                            ret[seq].Add(new BlastHashSeed((ushort)i, k));
                        }
                        else
                        {
                            List<BlastHashSeed> list = new List<BlastHashSeed>();
                            list.Add(new BlastHashSeed((ushort)i, k));
                            ret.Add(seq, list);
                        }

                    }

                }

                return RemoveLowComplexity(ret, nucCount);
            }

            public DNADictionary<List<BlastHashSeed>> GetBlastHashtable(NucleotideSequence databaseSequence)
            {
                DNADictionary<List<BlastHashSeed>> ret = new DNADictionary<List<BlastHashSeed>>(wordsize);
                ISequenceContainer template = SequenceFactory.GetDefaultSequence();


                ISequenceContainer data = databaseSequence.Sequence;
                for (int k = 0; k < data.Length - wordsize + 1; k++)
                {
                    ISequenceContainer seq = data.ConvertTo(template,k, wordsize);
                    if (ret.ContainsKey(seq))
                    {
                        ret[seq].Add(new BlastHashSeed(k));
                    }
                    else
                    {
                        List<BlastHashSeed> list = new List<BlastHashSeed>();
                        list.Add(new BlastHashSeed(k));
                        ret.Add(seq, list);
                    }

                }
                return RemoveLowComplexity(ret, data.Length);
            }

            /// <summary>
            /// Retrieve the wordsize
            /// </summary>
            public int Wordsize
            {
                get
                {
                    return this.wordsize;
                }
            }

            public int Stepsize
            {
                get
                {
                    return 1;
                }
            }

            /// <summary>
            /// Remove low complexity regions from the hash table,
            /// low complexity regions slow down the algorithm and create faulty alignments
            /// </summary>
            /// <param name="toRemove"></param>
            /// <param name="nucCount"></param>
            /// <returns></returns>
            private DNADictionary<List<BlastHashSeed>> RemoveLowComplexity(DNADictionary<List<BlastHashSeed>> toRemove, long nucCount)
            {
                int entries = (int)Math.Pow(4.0, this.wordsize);
                double expected = ((double)(nucCount - wordsize)) / (double)entries;

                int thresholdCount = (int)(expected * this.lowComplexityCutoff);
                if (thresholdCount < 5) thresholdCount = 5;

                for (int i = 0; i < entries; i++)
                {
                    List<BlastHashSeed> l = toRemove.GetValueForHash(i);
                    if (l != null)
                    {
                        if (toRemove.GetValueForHash(i).Count >= thresholdCount) toRemove[i] = null;
                        else l.TrimExcess(); //If the entry should be retained, trim excessively used capacity of the array!
                    }
                }

                return toRemove;

            }


        }


        /// <summary>
        /// Determines whether a given alignment is significant, soleley considering its score.
        /// The score of the alignment has to be greater or equal to the minimum score;
        /// </summary>
        public class BlastSignificator_Score : ISignificator
        {
            private int minScore;
            public BlastSignificator_Score(int minScore)
            {
                this.minScore = minScore;
            }







            public bool IsSignificant(IPairwiseAlignmentContainer simpleAlignment)
            {
                if (simpleAlignment.Score >= minScore) return true;
                else return false;
            }


        }

        /// <summary>
        /// 
        /// </summary>
        /// <summary>
        /// Determines whether a given alignment is significant, All are significant

        /// </summary>
        public class BlastSignificator_All : ISignificator
        {

            public BlastSignificator_All()
            {

            }







            public bool IsSignificant(IPairwiseAlignmentContainer simpleAlignment)
            {

                return true;
            }


        }


   
        /// <summary>
        /// Creates an alignment for a specified position in the database and in the query sequence.
        /// The alignment is extended in the 5' and in the 3' direction, in a stepwise manner and dynamic manner, ie. the alignment window is
        /// adjusted with gaps.
        /// </summary>
        public class AnchoredDynamicBandedSmithWatermanGotoh : IAnchoredDynamicProgramming
        {
            //Working variables
            private ISequenceContainer databaseSequence;
            private ISequenceContainer querySequence;
            private int smithWatermanSize;
            private int overlap;
            private SubstitutionMatrix substitutionMatrix;
            private PairwiseAlignment alignment;
            private int anchorPosDatabase;
            private int anchorPosQuery;
            private int activationThreshold;

            //Resulting variables
            private int alignmentStartDatabase = 0;
            private int alignmentStartQuery = 0;
            private int alignmentEndDatabase = 0;
            private int alignmentEndQuery = 0;
            private float highscore = 0.0F;

            //Internal control flags;
            private bool? usedDynmaicBanded = null;

            public AnchoredDynamicBandedSmithWatermanGotoh(int smithWatermanSize, int overhead, int activationThreshold)
            {
                this.smithWatermanSize = smithWatermanSize;
                this.overlap = overhead;
                this.activationThreshold = activationThreshold;
            }

            public AnchoredDynamicBandedSmithWatermanGotoh(ISequenceContainer databaseSequence, ISequenceContainer querySequence, SubstitutionMatrix substitutionMatrix, int anchorPositionDatabase, int anchorPositionQuery, int smithWatermanSize, int overhead, int activationThreshold)
            {
                this.databaseSequence = databaseSequence;
                this.querySequence = querySequence;
                this.anchorPosDatabase = anchorPositionDatabase - 1; //Convert into zero based array
                this.anchorPosQuery = anchorPositionQuery - 1;       //Convert into zero based array
                this.substitutionMatrix = substitutionMatrix;
                this.smithWatermanSize = smithWatermanSize;
                this.activationThreshold = activationThreshold;
                this.overlap = overhead;
                CreateAlignment();
            }


            private void CreateAlignment()
            {
                if (querySequence.Length > activationThreshold)
                {
                    CreateAlignmentDynamic();

                }
                else CreateAlignmentCommon();

                if (alignment.Length < 1) ResetValues();
            }

            private void ResetValues()
            {
                this.alignment = null;
                this.highscore = 0;
                this.alignmentEndDatabase = 0;
                this.alignmentEndQuery = 0;
                this.alignmentStartDatabase = 0;
                this.alignmentStartQuery = 0;
                this.highscore = 0;
            }

            private void CreateAlignmentCommon()
            {
                AnchoredSmithWatermanGotoh asmg = new AnchoredSmithWatermanGotoh(this.databaseSequence, this.querySequence, this.substitutionMatrix, this.anchorPosDatabase, this.anchorPosQuery);
                this.alignment = asmg.GetAlignment();
                this.alignmentEndDatabase = asmg.End_Database;
                this.alignmentEndQuery = asmg.End_Query;
                this.alignmentStartDatabase = asmg.Start_Database;
                this.alignmentStartQuery = asmg.Start_Query;
                this.highscore = asmg.Score;
                this.usedDynmaicBanded = false;

            }

            private void CreateAlignmentDynamic()
            {
                PairwiseAlignmentBuilder builder = new PairwiseAlignmentBuilder();
                ISequenceContainer template = SequenceFactory.GetDefaultSmithWatermanSequence();
                bool iterate = true;

                ///5 Prime EXTENSION
                ///Array:012345
                ///Seq:  AAAAAA eg: anchor at 5 -> smith waterman from position=0, length=5
                ///01234567
                ///AATTTAAT anchor at 7, smithWatermanLength=5 -> start at position 2, length 5
                int iterateDatabasePos = anchorPosDatabase;
                int iterateQueryPos = anchorPosQuery;
                while (iterate)
                {
                    //Get the sequence for which the Smith Waterman alignment should be constructed
                    ISequenceContainer databaseFragment;
                    ISequenceContainer queryFragment;
                    if (iterateDatabasePos - smithWatermanSize <= 0 && iterateQueryPos - smithWatermanSize > 0)
                    {
                        databaseFragment = databaseSequence.ConvertTo(template,0, iterateDatabasePos);
                        queryFragment = querySequence.ConvertTo(template,iterateQueryPos - smithWatermanSize, smithWatermanSize);
                        iterate = false;
                    }
                    else if (iterateQueryPos - smithWatermanSize <= 0 && iterateDatabasePos - smithWatermanSize > 0)
                    {
                        databaseFragment = databaseSequence.ConvertTo(template,iterateDatabasePos - smithWatermanSize, smithWatermanSize);
                        queryFragment = querySequence.ConvertTo(template,0, iterateQueryPos);
                        iterate = false;
                    }
                    else if (iterateDatabasePos - smithWatermanSize <= 0 || iterateQueryPos - smithWatermanSize <= 0)
                    {
                        queryFragment = querySequence.ConvertTo(template,0, iterateQueryPos);
                        databaseFragment = databaseSequence.ConvertTo(template,0, iterateDatabasePos);
                        iterate = false;
                    }
                    else
                    {
                        databaseFragment = databaseSequence.ConvertTo(template,iterateDatabasePos - smithWatermanSize, smithWatermanSize);
                        queryFragment = querySequence.ConvertTo(template,iterateQueryPos - smithWatermanSize, smithWatermanSize);
                    }


                    if (databaseFragment.Length == 0 || queryFragment.Length == 0) break;

                    SmithWatermanGotoh_DynamicBanded_5p sm5 = new SmithWatermanGotoh_DynamicBanded_5p(databaseFragment, queryFragment, substitutionMatrix, iterate == false ? 0 : overlap);
                    //    -54321
                    //     TTTTT
                    //01234567890
                    //AAAAATTTTTN anchor is at 10, alignment starts at position -5, next anchor should be at position 5  -> 
                    iterateDatabasePos += sm5.Start_Database;
                    iterateQueryPos += sm5.Start_Query;
                    //append the new alignment at the 5 prime end of the builder
                    builder.Append_5_prime(sm5.GetAlignment());

                    this.highscore += sm5.Score;

                    if (sm5.EndofDynamicExtension)
                    {
                        iterate = false;
                    }
                }//End of iteration Highest scoring end of the 5 prime alignment was found
                //Set Alignmentbegin
                this.alignmentStartDatabase = iterateDatabasePos + 1;
                this.alignmentStartQuery = iterateQueryPos + 1;


                //3'-Extension
                iterateDatabasePos = anchorPosDatabase - 1;
                iterateQueryPos = anchorPosQuery - 1;
                iterate = true;
                //0123456
                //AAAAATT length=7, anchor=1, smithWatermanLength=5
                while (iterate)
                {
                    //Get the sequence for which the Smith Waterman alignment should be constructed
                    ISequenceContainer databaseFragment;
                    ISequenceContainer queryFragment;
                    if (iterateDatabasePos + smithWatermanSize >= databaseSequence.Length && iterateQueryPos + smithWatermanSize < querySequence.Length)
                    {
                        databaseFragment = databaseSequence.ConvertTo(template,iterateDatabasePos + 1, databaseSequence.Length - iterateDatabasePos);
                        queryFragment = querySequence.ConvertTo(template,iterateQueryPos + 1, smithWatermanSize);
                        iterate = false;
                    }
                    else if (iterateQueryPos + smithWatermanSize >= querySequence.Length && iterateDatabasePos + smithWatermanSize < databaseSequence.Length)
                    {
                        databaseFragment = databaseSequence.ConvertTo(template,iterateDatabasePos + 1, smithWatermanSize);
                        queryFragment = querySequence.ConvertTo(template,iterateQueryPos + 1, querySequence.Length - iterateQueryPos);
                        iterate = false;
                    }
                    else if (iterateDatabasePos + smithWatermanSize >= databaseSequence.Length || iterateQueryPos + smithWatermanSize >= querySequence.Length)
                    {
                        queryFragment = querySequence.ConvertTo(template,iterateQueryPos + 1, querySequence.Length - iterateQueryPos);
                        databaseFragment = databaseSequence.ConvertTo(template,iterateDatabasePos + 1, databaseSequence.Length - iterateDatabasePos);
                        iterate = false;
                    }
                    else
                    {
                        databaseFragment = databaseSequence.ConvertTo(template,iterateDatabasePos + 1, smithWatermanSize);
                        queryFragment = querySequence.ConvertTo(template,iterateQueryPos + 1, smithWatermanSize);
                    }

                    if (databaseFragment.Length == 0 || queryFragment.Length == 0) break;

                    SmithWatermanGotoh_DynamicBanded_3p sm3 = new SmithWatermanGotoh_DynamicBanded_3p(databaseFragment, queryFragment, substitutionMatrix, iterate == false ? 0 : overlap);
                    //    -54321
                    //     TTTTT
                    //01234567890
                    //AAAAATTTTTN anchor is at 10, alignment starts at position -5, next anchor should be at position 3  -> 
                    iterateDatabasePos += sm3.End_Database;
                    iterateQueryPos += sm3.End_Query;
                    //append the new alignment at the 5 prime end of the builder
                    builder.Append_3_prime(sm3.GetAlignment());

                    this.highscore += sm3.Score;
                    if (sm3.EndofDynamicExtension)
                    {
                        iterate = false;
                    }
                }//End of iteration, highest scoring segment at the 3' end has been found

                this.alignmentEndDatabase = iterateDatabasePos + 1;
                this.alignmentEndQuery = iterateQueryPos + 1;

                this.alignment = builder.GetAlignment();

                this.usedDynmaicBanded = true;

            }


            public int DynamicProgrammingMatrixLength
            {
                get
                {
                    return this.smithWatermanSize;
                }
            }
            public int DynamicProgrammingMatrixOverlapLength
            {
                get
                {

                    return this.overlap;
                }
            }
            public int ActivationThreshold
            {
                get
                {
                    return this.activationThreshold;
                }
            }

            public bool UsedDynamicBanding
            {
                get
                {
                    return usedDynmaicBanded.Value;
                }
            }




            #region ITwoSeqAlignment Members

            public PairwiseAlignment GetAlignment()
            {

                return alignment;
            }

            public float Score
            {
                get
                {
                    return this.highscore;
                }
            }

            public int Start_Database
            {
                get
                {

                    return this.alignmentStartDatabase;
                }
            }

            public int Start_Query
            {
                get
                {
                    return this.alignmentStartQuery;
                }
            }

            public int End_Database
            {
                get
                {
                    return this.alignmentEndDatabase;
                }
            }

            public int End_Query
            {
                get
                {
                    return this.alignmentEndQuery;
                }
            }
            public SubstitutionMatrix SubstitutionMatrix
            {
                get
                {
                    return this.substitutionMatrix;
                }
            }

            #endregion


            /// <summary>
            /// Prototype pattern, get new instance of the AnchoredStepwiseSmithWatermanGotoh.
            /// Parameters specific for AnchoredStepwiseSmithWatermanGotoh are copied from this prototype instance.
            /// Anchored dynamic programming algorithm may be used in an exchangeable fashin within BLAST implementations.
            /// </summary>
            /// <param name="database">the database sequence</param>
            /// <param name="query">the query sequence</param>
            /// <param name="substitutionMatrix">the substitution matrix</param>
            /// <param name="anchorPositionDatabase">the position of the anchor (word, k-tuple) with respect to the database sequence</param>
            /// <param name="anchorPositionQuery">the position of the anchro (word,k-tuple) with respect to the query sequence</param>
            /// <returns>a new instance of the dynamic programming</returns>
            public IAnchoredDynamicProgramming NewDynamicProgramming(ISequenceContainer database, ISequenceContainer query, SubstitutionMatrix substitutionMatrix, int anchorPositionDatabase, int anchorPositionQuery)
            {
                return new AnchoredDynamicBandedSmithWatermanGotoh(database, query, substitutionMatrix, anchorPositionDatabase, anchorPositionQuery, this.DynamicProgrammingMatrixLength, this.DynamicProgrammingMatrixOverlapLength, this.ActivationThreshold);
            }






            /// <summary>
            /// This overloaded method implementation yields the same results as the method without plus/plus orientation.
            /// In the AnchoredStepwiseSmithWatermanGotoh it does not make a difference whether the plus/plus or plus/minus orientation is used
            /// </summary>
            /// <param name="database"></param>
            /// <param name="query"></param>
            /// <param name="substitutionMatrix"></param>
            /// <param name="anchorPositionDatabase"></param>
            /// <param name="anchorPositionQuery"></param>
            /// <param name="plusplusOrientation"></param>
            /// <returns></returns>
            public IAnchoredDynamicProgramming NewDynamicProgramming(ISequenceContainer database, ISequenceContainer query, SubstitutionMatrix substitutionMatrix, int anchorPositionDatabase, int anchorPositionQuery, bool plusplusOrientation)
            {
                return this.NewDynamicProgramming(database, query, substitutionMatrix, anchorPositionDatabase, anchorPositionQuery);
            }


        }


        /// <summary>
        /// Creates an alignment for a specified position in the database and in the query sequence.
        /// The alignment is extended in the 5' and in the 3' direction, in a stepwise manner and dynamic manner, ie. the alignment window is
        /// adjusted with gaps.
        /// </summary>
        public class Anchored454DynamicBandedSmithWatermanGotoh : IAnchoredDynamicProgramming
        {
            //Working variables
            private ISequenceContainer databaseSequence;
            private ISequenceContainer querySequence;
            private int smithWatermanSize;
            private int overlap;
            private SubstitutionMatrix substitutionMatrix;
            private PairwiseAlignment alignment;
            private int anchorPosDatabase;
            private int anchorPosQuery;
            private int activationThreshold;
            private float homopolymerTransgressionPenalty;

            //Resulting variables
            private int alignmentStartDatabase = 0;
            private int alignmentStartQuery = 0;
            private int alignmentEndDatabase = 0;
            private int alignmentEndQuery = 0;
            private float highscore = 0.0F;
            private bool plusPlusOrientation;
            private bool? usedDynmaicBanded = null;


            public Anchored454DynamicBandedSmithWatermanGotoh(int smithWatermanSize, int overhead, int activationThreshold, float homopolymerTransgressionPenalty)
            {
                this.smithWatermanSize = smithWatermanSize;
                this.overlap = overhead;
                this.activationThreshold = activationThreshold;
                this.homopolymerTransgressionPenalty = homopolymerTransgressionPenalty;
            }

            public Anchored454DynamicBandedSmithWatermanGotoh(ISequenceContainer databaseSequence, ISequenceContainer querySequence, SubstitutionMatrix substitutionMatrix,bool plusPlusOrientation, int anchorPositionDatabase, int anchorPositionQuery, int smithWatermanSize, int overhead, int activationThreshold,float homopolymerTransgressionPenalty)
            {
                this.databaseSequence = databaseSequence;
                this.querySequence = querySequence;
                this.anchorPosDatabase = anchorPositionDatabase - 1; //Convert into zero based array
                this.anchorPosQuery = anchorPositionQuery - 1;       //Convert into zero based array
                this.substitutionMatrix = substitutionMatrix;
                this.smithWatermanSize = smithWatermanSize;
                this.activationThreshold = activationThreshold;
                this.plusPlusOrientation = plusPlusOrientation;
                this.homopolymerTransgressionPenalty = homopolymerTransgressionPenalty;
                this.overlap = overhead;

                CreateAlignment();
            }


            private void CreateAlignment()
            {
                if (querySequence.Length > activationThreshold)
                {
                    if (plusPlusOrientation) CreateAlignmentDynamicPlus();
                    else CreateAlignmentDynamicMinus(); 
                }
                else CreateAlignmentCommon();

                if (alignment.Length < 1) ResetValues();
            }

            private void ResetValues()
            {
                this.alignment = null;
                this.highscore = 0;
                this.alignmentEndDatabase = 0;
                this.alignmentEndQuery = 0;
                this.alignmentStartDatabase = 0;
                this.alignmentStartQuery = 0;
                this.highscore = 0;
            }

            private void CreateAlignmentCommon()
            {
                Anchored454SmithWatermanGotoh asmg = new Anchored454SmithWatermanGotoh(this.databaseSequence, this.querySequence, this.substitutionMatrix, this.anchorPosDatabase, this.anchorPosQuery,plusPlusOrientation,homopolymerTransgressionPenalty);
                this.alignment = asmg.GetAlignment();
                this.alignmentEndDatabase = asmg.End_Database;
                this.alignmentEndQuery = asmg.End_Query;
                this.alignmentStartDatabase = asmg.Start_Database;
                this.alignmentStartQuery = asmg.Start_Query;
                this.highscore = asmg.Score;
                this.usedDynmaicBanded = false;
            }


            private void CreateAlignmentDynamicPlus()
            {
                PairwiseAlignmentBuilder builder = new PairwiseAlignmentBuilder();
                ISequenceContainer template = SequenceFactory.GetDefaultSmithWatermanSequence();
                bool iterate = true;

                ///5 Prime EXTENSION
                ///Array:012345
                ///Seq:  AAAAAA eg: anchor at 5 -> smith waterman from position=0, length=5
                ///01234567
                ///AATTTAAT anchor at 7, smithWatermanLength=5 -> start at position 2, length 5
                int iterateDatabasePos = anchorPosDatabase;
                int iterateQueryPos = anchorPosQuery;
                while (iterate)
                {
                    //Get the sequence for which the Smith Waterman alignment should be constructed
                    ISequenceContainer databaseFragment;
                    ISequenceContainer queryFragment;
                    if (iterateDatabasePos - smithWatermanSize <= 0 && iterateQueryPos - smithWatermanSize > 0)
                    {
                        databaseFragment = databaseSequence.ConvertTo(template, 0, iterateDatabasePos);
                        queryFragment = querySequence.ConvertTo(template, iterateQueryPos - smithWatermanSize, smithWatermanSize);
                        iterate = false;
                    }
                    else if (iterateQueryPos - smithWatermanSize <= 0 && iterateDatabasePos - smithWatermanSize > 0)
                    {
                        databaseFragment = databaseSequence.ConvertTo(template, iterateDatabasePos - smithWatermanSize, smithWatermanSize);
                        queryFragment = querySequence.ConvertTo(template, 0, iterateQueryPos);
                        iterate = false;
                    }
                    else if (iterateDatabasePos - smithWatermanSize <= 0 || iterateQueryPos - smithWatermanSize <= 0)
                    {
                        queryFragment = querySequence.ConvertTo(template, 0, iterateQueryPos);
                        databaseFragment = databaseSequence.ConvertTo(template, 0, iterateDatabasePos);
                        iterate = false;
                    }
                    else
                    {
                        databaseFragment = databaseSequence.ConvertTo(template, iterateDatabasePos - smithWatermanSize, smithWatermanSize);
                        queryFragment = querySequence.ConvertTo(template, iterateQueryPos - smithWatermanSize, smithWatermanSize);
                    }


                    if (databaseFragment.Length == 0 || queryFragment.Length == 0) break;

                    SmithWatermanGotoh_DynamicBanded_454P_5p sm5 = new SmithWatermanGotoh_DynamicBanded_454P_5p(databaseFragment, queryFragment, substitutionMatrix,homopolymerTransgressionPenalty ,iterate == false ? 0 : overlap);
                    
                    //    -54321
                    //     TTTTT
                    //01234567890
                    //AAAAATTTTTN anchor is at 10, alignment starts at position -5, next anchor should be at position 5  -> 
                    iterateDatabasePos += sm5.Start_Database;
                    iterateQueryPos += sm5.Start_Query;
                    //append the new alignment at the 5 prime end of the builder
                    builder.Append_5_prime(sm5.GetAlignment());

                    this.highscore += sm5.Score;

                    if (sm5.EndofDynamicExtension)
                    {
                        iterate = false;
                    }
                }//End of iteration Highest scoring end of the 5 prime alignment was found
                //Set Alignmentbegin
                this.alignmentStartDatabase = iterateDatabasePos + 1;
                this.alignmentStartQuery = iterateQueryPos + 1;


                //3'-Extension
                iterateDatabasePos = anchorPosDatabase - 1;
                iterateQueryPos = anchorPosQuery - 1;
                iterate = true;
                //0123456
                //AAAAATT length=7, anchor=1, smithWatermanLength=5
                while (iterate)
                {
                    //Get the sequence for which the Smith Waterman alignment should be constructed
                    ISequenceContainer databaseFragment;
                    ISequenceContainer queryFragment;
                    if (iterateDatabasePos + smithWatermanSize >= databaseSequence.Length && iterateQueryPos + smithWatermanSize < querySequence.Length)
                    {
                        databaseFragment = databaseSequence.ConvertTo(template, iterateDatabasePos + 1, databaseSequence.Length - iterateDatabasePos);
                        queryFragment = querySequence.ConvertTo(template, iterateQueryPos + 1, smithWatermanSize);
                        iterate = false;
                    }
                    else if (iterateQueryPos + smithWatermanSize >= querySequence.Length && iterateDatabasePos + smithWatermanSize < databaseSequence.Length)
                    {
                        databaseFragment = databaseSequence.ConvertTo(template, iterateDatabasePos + 1, smithWatermanSize);
                        queryFragment = querySequence.ConvertTo(template, iterateQueryPos + 1, querySequence.Length - iterateQueryPos);
                        iterate = false;
                    }
                    else if (iterateDatabasePos + smithWatermanSize >= databaseSequence.Length || iterateQueryPos + smithWatermanSize >= querySequence.Length)
                    {
                        queryFragment = querySequence.ConvertTo(template, iterateQueryPos + 1, querySequence.Length - iterateQueryPos);
                        databaseFragment = databaseSequence.ConvertTo(template, iterateDatabasePos + 1, databaseSequence.Length - iterateDatabasePos);
                        iterate = false;
                    }
                    else
                    {
                        databaseFragment = databaseSequence.ConvertTo(template, iterateDatabasePos + 1, smithWatermanSize);
                        queryFragment = querySequence.ConvertTo(template, iterateQueryPos + 1, smithWatermanSize);
                    }

                    if (databaseFragment.Length == 0 || queryFragment.Length == 0) break;

                    SmithWatermanGotoh_DynamicBanded_454P_3p sm3 = new SmithWatermanGotoh_DynamicBanded_454P_3p(databaseFragment, queryFragment, substitutionMatrix, iterate == false ? 0 : overlap);
                    //    -54321
                    //     TTTTT
                    //01234567890
                    //AAAAATTTTTN anchor is at 10, alignment starts at position -5, next anchor should be at position 3  -> 
                    iterateDatabasePos += sm3.End_Database;
                    iterateQueryPos += sm3.End_Query;
                    //append the new alignment at the 5 prime end of the builder
                    builder.Append_3_prime(sm3.GetAlignment());

                    this.highscore += sm3.Score;
                    if (sm3.EndofDynamicExtension)
                    {
                        iterate = false;
                    }
                }//End of iteration, highest scoring segment at the 3' end has been found

                this.alignmentEndDatabase = iterateDatabasePos + 1;
                this.alignmentEndQuery = iterateQueryPos + 1;

                this.alignment = builder.GetAlignment();
                usedDynmaicBanded = true;
            }


            private void CreateAlignmentDynamicMinus()
            {
                PairwiseAlignmentBuilder builder = new PairwiseAlignmentBuilder();
                ISequenceContainer template = SequenceFactory.GetDefaultSmithWatermanSequence();
                bool iterate = true;

                ///5 Prime EXTENSION
                ///Array:012345
                ///Seq:  AAAAAA eg: anchor at 5 -> smith waterman from position=0, length=5
                ///01234567
                ///AATTTAAT anchor at 7, smithWatermanLength=5 -> start at position 2, length 5
                int iterateDatabasePos = anchorPosDatabase;
                int iterateQueryPos = anchorPosQuery;
                while (iterate)
                {
                    //Get the sequence for which the Smith Waterman alignment should be constructed
                    ISequenceContainer databaseFragment;
                    ISequenceContainer queryFragment;
                    if (iterateDatabasePos - smithWatermanSize <= 0 && iterateQueryPos - smithWatermanSize > 0)
                    {
                        databaseFragment = SequenceUtility.ReverseSequence(databaseSequence.ConvertTo(template, 0, iterateDatabasePos));
                        queryFragment =  SequenceUtility.ReverseSequence(querySequence.ConvertTo(template, iterateQueryPos - smithWatermanSize, smithWatermanSize));
                        iterate = false;
                    }
                    else if (iterateQueryPos - smithWatermanSize <= 0 && iterateDatabasePos - smithWatermanSize > 0)
                    {
                        databaseFragment = SequenceUtility.ReverseSequence(databaseSequence.ConvertTo(template, iterateDatabasePos - smithWatermanSize, smithWatermanSize));
                        queryFragment = SequenceUtility.ReverseSequence(querySequence.ConvertTo(template, 0, iterateQueryPos));
                        iterate = false;
                    }
                    else if (iterateDatabasePos - smithWatermanSize <= 0 || iterateQueryPos - smithWatermanSize <= 0)
                    {
                        queryFragment = SequenceUtility.ReverseSequence(querySequence.ConvertTo(template, 0, iterateQueryPos));
                        databaseFragment = SequenceUtility.ReverseSequence(databaseSequence.ConvertTo(template, 0, iterateDatabasePos));
                        iterate = false;
                    }
                    else
                    {
                        databaseFragment = SequenceUtility.ReverseSequence(databaseSequence.ConvertTo(template, iterateDatabasePos - smithWatermanSize, smithWatermanSize));
                        queryFragment = SequenceUtility.ReverseSequence(querySequence.ConvertTo(template, iterateQueryPos - smithWatermanSize, smithWatermanSize));
                    }


                    if (databaseFragment.Length == 0 || queryFragment.Length == 0) break;
                    //Carefull this is the 5 prime extension but since we are here using here reversed sequences the 3p extending SmithWaterman is used
                    SmithWatermanGotoh_DynamicBanded_454P_3p sm5 = new SmithWatermanGotoh_DynamicBanded_454P_3p(databaseFragment, queryFragment, substitutionMatrix,this.homopolymerTransgressionPenalty, iterate == false ? 0 : overlap);
                    //    -54321
                    //     TTTTT
                    //01234567890
                    //AAAAATTTTTN anchor is at 10, alignment starts at position -5, next anchor should be at position 5  -> 
                    //Again be careful; since the sequences have been reverssed
                    iterateDatabasePos += -sm5.End_Database;
                    iterateQueryPos += -sm5.End_Query;
                    builder.Append_5_prime(PairwiseAlignmentUtility.ReverseAlignment(sm5.GetAlignment()));

                    this.highscore += sm5.Score;

                    if (sm5.EndofDynamicExtension)
                    {
                        iterate = false;
                    }
                }//End of iteration Highest scoring end of the 5 prime alignment was found
                //Set Alignmentbegin
                this.alignmentStartDatabase = iterateDatabasePos + 1;
                this.alignmentStartQuery = iterateQueryPos + 1;


                //3'-Extension
                iterateDatabasePos = anchorPosDatabase - 1;
                iterateQueryPos = anchorPosQuery - 1;
                iterate = true;
                //0123456
                //AAAAATT length=7, anchor=1, smithWatermanLength=5
                while (iterate)
                {
                    //Get the sequence for which the Smith Waterman alignment should be constructed
                    ISequenceContainer databaseFragment;
                    ISequenceContainer queryFragment;
                    if (iterateDatabasePos + smithWatermanSize >= databaseSequence.Length && iterateQueryPos + smithWatermanSize < querySequence.Length)
                    {
                        databaseFragment = SequenceUtility.ReverseSequence(databaseSequence.ConvertTo(template, iterateDatabasePos + 1, databaseSequence.Length - iterateDatabasePos));
                        queryFragment = SequenceUtility.ReverseSequence(querySequence.ConvertTo(template, iterateQueryPos + 1, smithWatermanSize));
                        iterate = false;
                    }
                    else if (iterateQueryPos + smithWatermanSize >= querySequence.Length && iterateDatabasePos + smithWatermanSize < databaseSequence.Length)
                    {
                        databaseFragment = SequenceUtility.ReverseSequence(databaseSequence.ConvertTo(template, iterateDatabasePos + 1, smithWatermanSize));
                        queryFragment = SequenceUtility.ReverseSequence(querySequence.ConvertTo(template, iterateQueryPos + 1, querySequence.Length - iterateQueryPos));
                        iterate = false;
                    }
                    else if (iterateDatabasePos + smithWatermanSize >= databaseSequence.Length || iterateQueryPos + smithWatermanSize >= querySequence.Length)
                    {
                        queryFragment = SequenceUtility.ReverseSequence(querySequence.ConvertTo(template, iterateQueryPos + 1, querySequence.Length - iterateQueryPos));
                        databaseFragment = SequenceUtility.ReverseSequence(databaseSequence.ConvertTo(template, iterateDatabasePos + 1, databaseSequence.Length - iterateDatabasePos));
                        iterate = false;
                    }
                    else
                    {
                        databaseFragment = SequenceUtility.ReverseSequence(databaseSequence.ConvertTo(template, iterateDatabasePos + 1, smithWatermanSize));
                        queryFragment = SequenceUtility.ReverseSequence(querySequence.ConvertTo(template, iterateQueryPos + 1, smithWatermanSize));
                    }

                    if (databaseFragment.Length == 0 || queryFragment.Length == 0) break;

                    SmithWatermanGotoh_DynamicBanded_454P_5p sm3 = new SmithWatermanGotoh_DynamicBanded_454P_5p(databaseFragment, queryFragment, substitutionMatrix,this.homopolymerTransgressionPenalty, iterate == false ? 0 : overlap);

                    //    -54321
                    //     TTTTT
                    //01234567890
                    //AAAAATTTTTN anchor is at 10, alignment starts at position -5, next anchor should be at position 3  -> 
                    iterateDatabasePos += -sm3.Start_Database;
                    iterateQueryPos += -sm3.Start_Query;
                    //append the new alignment at the 5 prime end of the builder
                    builder.Append_3_prime(PairwiseAlignmentUtility.ReverseAlignment(sm3.GetAlignment()));

                    this.highscore += sm3.Score;
                    if (sm3.EndofDynamicExtension)
                    {
                        iterate = false;
                    }
                }//End of iteration, highest scoring segment at the 3' end has been found

                this.alignmentEndDatabase = iterateDatabasePos + 1;
                this.alignmentEndQuery = iterateQueryPos + 1;

                this.alignment = builder.GetAlignment();

                this.usedDynmaicBanded = true;

            }



            public int DynamicProgrammingMatrixLength
            {
                get
                {
                    return this.smithWatermanSize;
                }
            }
            public int DynamicProgrammingMatrixOverlapLength
            {
                get
                {

                    return this.overlap;
                }
            }
            public int ActivationThreshold
            {
                get
                {
                    return this.activationThreshold;
                }
            }

            public bool UsedDynamicBanding
            {
                get
                {
                    return usedDynmaicBanded.Value;
                }
            }




            #region ITwoSeqAlignment Members

            public PairwiseAlignment GetAlignment()
            {
                return alignment;
            }

            public float Score
            {
                get
                {
                    return this.highscore;
                }
            }

            public int Start_Database
            {
                get
                {

                    return this.alignmentStartDatabase;
                }
            }

            public int Start_Query
            {
                get
                {
                    return this.alignmentStartQuery;
                }
            }

            public int End_Database
            {
                get
                {
                    return this.alignmentEndDatabase;
                }
            }

            public int End_Query
            {
                get
                {
                    return this.alignmentEndQuery;
                }
            }
            public SubstitutionMatrix SubstitutionMatrix
            {
                get
                {
                    return this.substitutionMatrix;
                }
            }

            #endregion


            /// <summary>
            /// Prototype pattern, get new instance of the AnchoredStepwiseSmithWatermanGotoh.
            /// Parameters specific for AnchoredStepwiseSmithWatermanGotoh are copied from this prototype instance.
            /// Anchored dynamic programming algorithm may be used in an exchangeable fashin within BLAST implementations.
            /// </summary>
            /// <param name="database">the database sequence</param>
            /// <param name="query">the query sequence</param>
            /// <param name="substitutionMatrix">the substitution matrix</param>
            /// <param name="anchorPositionDatabase">the position of the anchor (word, k-tuple) with respect to the database sequence</param>
            /// <param name="anchorPositionQuery">the position of the anchro (word,k-tuple) with respect to the query sequence</param>
            /// <returns>a new instance of the dynamic programming</returns>
            public IAnchoredDynamicProgramming NewDynamicProgramming(ISequenceContainer database, ISequenceContainer query, SubstitutionMatrix substitutionMatrix, int anchorPositionDatabase, int anchorPositionQuery)
            {
                return new AnchoredDynamicBandedSmithWatermanGotoh(database, query, substitutionMatrix, anchorPositionDatabase, anchorPositionQuery, this.DynamicProgrammingMatrixLength, this.DynamicProgrammingMatrixOverlapLength, this.ActivationThreshold);
            }






            /// <summary>
            /// This overloaded method implementation yields the same results as the method without plus/plus orientation.
            /// In the AnchoredStepwiseSmithWatermanGotoh it does not make a difference whether the plus/plus or plus/minus orientation is used
            /// </summary>
            /// <param name="database"></param>
            /// <param name="query"></param>
            /// <param name="substitutionMatrix"></param>
            /// <param name="anchorPositionDatabase"></param>
            /// <param name="anchorPositionQuery"></param>
            /// <param name="plusplusOrientation"></param>
            /// <returns></returns>
            public IAnchoredDynamicProgramming NewDynamicProgramming(ISequenceContainer database, ISequenceContainer query, SubstitutionMatrix substitutionMatrix, int anchorPositionDatabase, int anchorPositionQuery, bool plusplusOrientation)
            {
                return this.NewDynamicProgramming(database, query, substitutionMatrix, anchorPositionDatabase, anchorPositionQuery);
            }


        }





        /// <summary>
        /// Creates an ungapped alignment, 
        /// </summary>
        public class AnchoredUngappedExtension : IAnchoredDynamicProgramming
        {
            //Working variables
            private ISequenceContainer databaseSequence;
            private ISequenceContainer querySequence;
            private SubstitutionMatrix substitutionMatrix;
            private PairwiseAlignment alignment;
            private int anchorPosDatabase;
            private int anchorPosQuery;
            private int testLength;

            //Resulting variables
            private int alignmentStartDatabase = 0;
            private int alignmentStartQuery = 0;
            private int alignmentEndDatabase = 0;
            private int alignmentEndQuery = 0;
            private float highscore = 0.0F;

            public AnchoredUngappedExtension(ISequenceContainer databaseSequence, ISequenceContainer querySequence, SubstitutionMatrix substitutionMatrix, int anchorPositionDatabase, int anchorPositionQuery, int testLength)
            {
                this.databaseSequence = databaseSequence;
                this.querySequence = querySequence;
                this.anchorPosDatabase = anchorPositionDatabase - 1; //Convert into zero based array
                this.anchorPosQuery = anchorPositionQuery - 1;       //Convert into zero based array
                this.substitutionMatrix = substitutionMatrix;

                this.testLength = testLength;

                CreateAlignment();
            }

            private void CreateAlignment()
            {


                int iterateDatabase = anchorPosDatabase;
                int iterateQuery = anchorPosQuery;
                int separation = 0;
                float score = highscore;

                //left extension
                while (iterateDatabase >= 0 && iterateQuery >= 0 && separation <= testLength)
                {
                    score += substitutionMatrix[databaseSequence[iterateDatabase], querySequence[iterateQuery]];
                    if (score > highscore)
                    {
                        highscore = score;
                        separation = 0;
                        this.alignmentStartDatabase = iterateDatabase + 1;
                        this.alignmentStartQuery = iterateQuery + 1;
                    }
                    iterateDatabase--;
                    iterateQuery--;
                    separation++;
                }

                //Right extension
                iterateDatabase = anchorPosDatabase + 1;
                iterateQuery = anchorPosQuery + 1;
                score = highscore;
                while (iterateDatabase < databaseSequence.Length && iterateQuery < querySequence.Length && separation <= testLength)
                {
                    score += substitutionMatrix[databaseSequence[iterateDatabase], querySequence[iterateQuery]];
                    if (score > highscore)
                    {
                        highscore = score;
                        separation = 0;
                        this.alignmentEndDatabase = iterateDatabase + 1;
                        this.alignmentStartQuery = iterateQuery + 1;
                    }
                    iterateDatabase++;
                    iterateQuery++;
                    separation++;
                }




            }







            #region ITwoSeqAlignment Members

            public PairwiseAlignment GetAlignment()
            {

                if (alignment == null)
                {
                    this.alignment = new PairwiseAlignment(databaseSequence.SubSequence(alignmentStartDatabase - 1, alignmentEndDatabase - 1), querySequence.SubSequence(alignmentStartQuery, alignmentEndQuery));
                }
                return alignment;
            }

            public float Score
            {
                get
                {
                    return this.highscore;
                }
            }

            public int Start_Database
            {
                get
                {

                    return this.alignmentStartDatabase;
                }
            }

            public int Start_Query
            {
                get
                {
                    return this.alignmentStartQuery;
                }
            }

            public int End_Database
            {
                get
                {
                    return this.alignmentEndDatabase;
                }
            }

            public int End_Query
            {
                get
                {
                    return this.alignmentEndQuery;
                }
            }
            public SubstitutionMatrix SubstitutionMatrix
            {
                get
                {
                    return this.substitutionMatrix;
                }
            }

            #endregion

            public int LengthOfUngapedSimilarityTest
            {
                get
                {
                    return this.testLength;
                }
            }
            /// <summary>
            /// Prototype pattern, creates a new instance of AnchoredUngapedExtension using the specified parameters.
            /// Parameters specific for AnchoredUngapedExtension are copied from this prototype instance.
            /// This implementation actually does not use dynamic programming algorithm, it just extends an anchor in both directions.
            /// May be used in BLAST implementations.
            /// </summary>
            /// <param name="database"></param>
            /// <param name="query"></param>
            /// <param name="substitutionMatrix"></param>
            /// <param name="anchorPositionDatabase"></param>
            /// <param name="anchorPositionQuery"></param>
            /// <returns></returns>
            public IAnchoredDynamicProgramming NewDynamicProgramming(ISequenceContainer database, ISequenceContainer query, SubstitutionMatrix substitutionMatrix, int anchorPositionDatabase, int anchorPositionQuery)
            {
                return new AnchoredUngappedExtension(database, query, substitutionMatrix, anchorPositionDatabase, anchorPositionQuery, this.LengthOfUngapedSimilarityTest);
            }




            /// <summary>
            /// This overloaded method implementation yields the same results as the version without plus/plus orientation. 
            /// With the AnchoredUngappedExtension there is no difference whether the plus/plus or the plus/minus orientation is used
            /// </summary>
            /// <param name="database"></param>
            /// <param name="query"></param>
            /// <param name="substitutionMatrix"></param>
            /// <param name="anchorPositionDatabase"></param>
            /// <param name="anchorPositionQuery"></param>
            /// <param name="plusplusOrientation"></param>
            /// <returns></returns>
            public IAnchoredDynamicProgramming NewDynamicProgramming(ISequenceContainer database, ISequenceContainer query, SubstitutionMatrix substitutionMatrix, int anchorPositionDatabase, int anchorPositionQuery, bool plusplusOrientation)
            {
                return this.NewDynamicProgramming(database, query, substitutionMatrix, anchorPositionDatabase, anchorPositionQuery);
            }


        }

        /// <summary>
        /// Creates a SmithWatermanGotoh alignment for the specified region
        /// </summary>
        public class Anchored454SmithWatermanGotoh : IAnchoredDynamicProgramming
        {
            //Working variables
            private ISequenceContainer databaseSequence;
            private ISequenceContainer querySequence;
            private SubstitutionMatrix substitutionMatrix;
            private PairwiseAlignment alignment;
            private int anchorPosDatabase;
            private int anchorPosQuery;

            private float homopolymereTransgressionBoundary;

            //boundary variables
            private int? maxBoundary_5P = 50;
            private int? maxBoundary_3P = 200;

            //Resulting variables
            private int alignmentStartDatabase = 0;
            private int alignmentStartQuery = 0;
            private int alignmentEndDatabase = 0;
            private int alignmentEndQuery = 0;
            private float highscore = 0.0F;
            private bool plusplusOrientation = true;
            private int subDatabaseStart = 0;
            

            public Anchored454SmithWatermanGotoh(ISequenceContainer databaseSequence, ISequenceContainer querySequence, SubstitutionMatrix substitutionMatrix, int anchorPositionDatabase, int anchorPositionQuery, bool plusplusOrientation, float homopolymereTransgressionBoundary)
            {
                this.databaseSequence = databaseSequence;
                this.querySequence = querySequence;
                this.substitutionMatrix = substitutionMatrix;
                this.anchorPosDatabase = anchorPositionDatabase;
                this.anchorPosQuery = anchorPositionQuery;
                this.plusplusOrientation = plusplusOrientation;
                this.homopolymereTransgressionBoundary = homopolymereTransgressionBoundary;

                InitiateAlignment();

            }

            /// <summary>
            /// This constructor should be used if the implementation is going to serve as prototype
            /// </summary>
            /// <param name="overheadInPercent"></param>
            public Anchored454SmithWatermanGotoh(float homopolymereTransgressionBoundary)
            {
                this.homopolymereTransgressionBoundary = homopolymereTransgressionBoundary;
            }

            public void InitiateAlignment()
            {
                ISequenceContainer template = SequenceFactory.GetDefaultSequence();
                

                //Calculate the overhead and restrict it when necessary to the limits
                int shift_5P = (int)((((double)anchorPosQuery * substitutionMatrix.HighestScore) - (double)substitutionMatrix.GapExistPenalty) / (double)substitutionMatrix.GapExtendPenalty);
                if (shift_5P < 10) shift_5P = 10;
                else if (maxBoundary_5P != null && shift_5P > maxBoundary_5P) shift_5P = maxBoundary_5P.Value;

                int shift_3P = (int)((((((double)querySequence.Length) / 2.0) * substitutionMatrix.HighestScore) - (double)substitutionMatrix.GapExistPenalty) / (double)substitutionMatrix.GapExtendPenalty);
                if (shift_3P < 10) shift_3P = 10;
                else if (maxBoundary_3P != null && shift_3P > maxBoundary_3P) shift_3P = maxBoundary_3P.Value;


                subDatabaseStart = (int)(anchorPosDatabase - anchorPosQuery - shift_5P); //Start of the subalignment; that is the small part in which teh sequences will be aligned
                if (subDatabaseStart < 0) subDatabaseStart = 0;
                int endOfDatabase = (int)(anchorPosDatabase + (querySequence.Length - anchorPosQuery) + shift_3P - 1.0F);
                if (endOfDatabase > databaseSequence.Length) endOfDatabase = databaseSequence.Length - 1;


                ISequenceContainer dataSub;
                dataSub = databaseSequence.ConvertTo(template, subDatabaseStart, endOfDatabase - subDatabaseStart + 1);


                IDynamicProgramming sw;
                if (plusplusOrientation) sw = new SmithWatermanGotoh_454P(dataSub, querySequence, substitutionMatrix, homopolymereTransgressionBoundary);
                else sw = new SmithWatermanGotoh_454M(dataSub, querySequence, substitutionMatrix, homopolymereTransgressionBoundary);

                this.alignment = sw.GetAlignment();
                this.alignmentStartDatabase = subDatabaseStart + sw.Start_Database;
                this.alignmentEndDatabase = subDatabaseStart + sw.End_Database;
                this.alignmentStartQuery = sw.Start_Query;
                this.alignmentEndQuery = sw.End_Query;
                this.highscore = sw.Score;
            }



            /// <summary>
            /// Prototype pattern, return new instance of Anchored454SmithWatermanGotoh
            /// </summary>
            /// <param name="database"></param>
            /// <param name="query"></param>
            /// <param name="substitutionMatrix"></param>
            /// <param name="anchorPositionDatabase"></param>
            /// <param name="anchorPositionQuery"></param>
            /// <returns></returns>
            public IAnchoredDynamicProgramming NewDynamicProgramming(ISequenceContainer database, ISequenceContainer query, SubstitutionMatrix substitutionMatrix, int anchorPositionDatabase, int anchorPositionQuery)
            {
                return new Anchored454SmithWatermanGotoh(database, query, substitutionMatrix, anchorPositionDatabase, anchorPositionQuery, this.plusplusOrientation, this.homopolymereTransgressionBoundary);
            }

            /// <summary>
            /// Return the alignment
            /// </summary>
            /// <returns></returns>
            public PairwiseAlignment GetAlignment()
            {
                return this.alignment;
            }

            public float Score
            {
                get
                {
                    return this.highscore;
                }
            }

            public SubstitutionMatrix SubstitutionMatrix
            {
                get
                {
                    return this.substitutionMatrix;
                }
            }

            public int Start_Database
            {
                get
                {
                    return this.alignmentStartDatabase;
                }
            }

            public int Start_Query
            {
                get
                {
                    return this.alignmentStartQuery;
                }
            }

            public int End_Database
            {
                get
                {
                    return this.alignmentEndDatabase;
                }
            }

            public int End_Query
            {
                get
                {
                    return this.alignmentEndQuery;
                }
            }


            public IAnchoredDynamicProgramming NewDynamicProgramming(ISequenceContainer database, ISequenceContainer query, SubstitutionMatrix substitutionMatrix, int anchorPositionDatabase, int anchorPositionQuery, bool plusplusOrientation)
            {
                return new Anchored454SmithWatermanGotoh(database, query, substitutionMatrix, anchorPositionDatabase, anchorPositionQuery, plusplusOrientation, this.homopolymereTransgressionBoundary);
            }


        }


        /// <summary>
        /// Creates a SmithWatermanGotoh alignment for the specified region
        /// </summary>
        public class AnchoredSmithWatermanGotoh : IAnchoredDynamicProgramming
        {
            //Working variables
            private ISequenceContainer databaseSequence;
            private ISequenceContainer querySequence;
            private SubstitutionMatrix substitutionMatrix;
            private PairwiseAlignment alignment;
            private int anchorPosDatabase;
            private int anchorPosQuery;

            //Changing variable;
            private int? maxBoundary_5P;
            private int? maxBoundary_3P;

            //Resulting variables
            private int alignmentStartDatabase = 0;
            private int alignmentStartQuery = 0;
            private int alignmentEndDatabase = 0;
            private int alignmentEndQuery = 0;
            private float highscore = 0.0F;
            private int subDatabaseStart = 0;

            public AnchoredSmithWatermanGotoh(ISequenceContainer databaseSequence, ISequenceContainer querySequence, SubstitutionMatrix substitutionMatrix, int anchorPositionDatabase, int anchorPositionQuery)
            {
                this.databaseSequence = databaseSequence;
                this.querySequence = querySequence;
                this.substitutionMatrix = substitutionMatrix;
                this.anchorPosDatabase = anchorPositionDatabase;
                this.anchorPosQuery = anchorPositionQuery;


                InitiateAlignment();

            }

            /// <summary>
            /// This constructor should be used if the implementation is going to serve as prototype
            /// </summary>
            /// <param name="overheadInPercent"></param>
            public AnchoredSmithWatermanGotoh()
            {

            }

            public void InitiateAlignment()
            {
                ISequenceContainer template = SequenceFactory.GetDefaultSequence();
                int shift_5P = (int)((((double)anchorPosQuery * substitutionMatrix.HighestScore) - (double)substitutionMatrix.GapExistPenalty) / (double)substitutionMatrix.GapExtendPenalty);
                if (shift_5P < 0) shift_5P = 10;
                else if (maxBoundary_5P != null && shift_5P > maxBoundary_5P) shift_5P = maxBoundary_5P.Value;

                int shift_3P = (int)((((((double)querySequence.Length) / 2.0) * substitutionMatrix.HighestScore) - (double)substitutionMatrix.GapExistPenalty) / (double)substitutionMatrix.GapExtendPenalty);
                if (shift_3P < 0) shift_3P = 10;
                else if (maxBoundary_3P != null && shift_3P > maxBoundary_3P) shift_3P = maxBoundary_3P.Value;

                subDatabaseStart = (int)(anchorPosDatabase - anchorPosQuery - shift_5P); //Start of the subalignment; that is the small part in which teh sequences will be aligned
                if (subDatabaseStart < 0) subDatabaseStart = 0;
                int endOfDatabase = (int)(anchorPosDatabase + (querySequence.Length - anchorPosQuery) + shift_3P - 1.0F);
                if (endOfDatabase > databaseSequence.Length) endOfDatabase = databaseSequence.Length - 1;

                ISequenceContainer dataSub;
                dataSub = databaseSequence.ConvertTo(template,subDatabaseStart, endOfDatabase - subDatabaseStart + 1);

                IDynamicProgramming sw = new SmithWatermanGotoh(dataSub, querySequence, substitutionMatrix);
                this.alignment = sw.GetAlignment();
                this.alignmentStartDatabase = subDatabaseStart + sw.Start_Database;
                this.alignmentEndDatabase = subDatabaseStart + sw.End_Database;
                this.alignmentStartQuery = sw.Start_Query;
                this.alignmentEndQuery = sw.End_Query;
                this.highscore = sw.Score;
            }



            /// <summary>
            /// Prototype pattern, return new instance of Anchored454SmithWatermanGotoh
            /// </summary>
            /// <param name="database"></param>
            /// <param name="query"></param>
            /// <param name="substitutionMatrix"></param>
            /// <param name="anchorPositionDatabase"></param>
            /// <param name="anchorPositionQuery"></param>
            /// <returns></returns>
            public IAnchoredDynamicProgramming NewDynamicProgramming(ISequenceContainer database, ISequenceContainer query, SubstitutionMatrix substitutionMatrix, int anchorPositionDatabase, int anchorPositionQuery)
            {
                return new AnchoredSmithWatermanGotoh(database, query, substitutionMatrix, anchorPositionDatabase, anchorPositionQuery);
            }

            /// <summary>
            /// Return the alignment
            /// </summary>
            /// <returns></returns>
            public PairwiseAlignment GetAlignment()
            {
                return this.alignment;
            }

            public float Score
            {
                get
                {
                    return this.highscore;
                }
            }

            public SubstitutionMatrix SubstitutionMatrix
            {
                get
                {
                    return this.substitutionMatrix;
                }
            }

            public int Start_Database
            {
                get
                {
                    return this.alignmentStartDatabase;
                }
            }

            public int Start_Query
            {
                get
                {
                    return this.alignmentStartQuery;
                }
            }

            public int End_Database
            {
                get
                {
                    return this.alignmentEndDatabase;
                }
            }

            public int End_Query
            {
                get
                {
                    return this.alignmentEndQuery;
                }
            }


            public IAnchoredDynamicProgramming NewDynamicProgramming(ISequenceContainer database, ISequenceContainer query, SubstitutionMatrix substitutionMatrix, int anchorPositionDatabase, int anchorPositionQuery, bool plusplusOrientation)
            {
                return new AnchoredSmithWatermanGotoh(database, query, substitutionMatrix, anchorPositionDatabase, anchorPositionQuery);
            }


        }



        /// <summary>
        /// Classes responsible for testing whether or not the specified alignment is significant have to implement this interface
        /// </summary>
        public interface ISignificator
        {
            bool IsSignificant(IPairwiseAlignmentContainer simpleAlignment);
        }

        public interface IBlastCore
        {
            List<IPairwiseAlignmentContainer> GetAlignments(List<NucleotideSequence> sequences);
        }

        /// <summary>
        /// Obtain a number of real seeds for a given List of BlastHashSeeds;
        /// BlastHashSeeds are in contrast to the BlastSeeds extremely memory efficient
        /// </summary>
        public interface IBlastSeedProcessor
        {
            List<BlastSeed> GetSeeds(List<BlastSeed> seeds);
        }

        public interface IBlastHashtableConstructor
        {
            DNADictionary<List<BlastHashSeed>> GetBlastHashtable(List<NucleotideSequence> databaseSequences);
            DNADictionary<List<BlastHashSeed>> GetBlastHashtable(NucleotideSequence databaseSequence);
            int Wordsize { get;}
            int Stepsize { get;}
        }


        public interface IBlastPostProcessor
        {
            List<IPairwiseAlignmentContainer> GetAlignments(List<PairwiseNucleotideSequenceAlignment> alignments);
        }

        /// <summary>
        /// Blast postprocessor; processes the blast search results
        /// Returns all hits unmodified
        /// </summary>
        public class BlastPostProcessorDefault : IBlastPostProcessor
        {
            public BlastPostProcessorDefault()
            {
            }
            public List<IPairwiseAlignmentContainer> GetAlignments(List<PairwiseNucleotideSequenceAlignment> alignments)
            {
                return alignments.ConvertAll<IPairwiseAlignmentContainer>(new Converter<PairwiseNucleotideSequenceAlignment, IPairwiseAlignmentContainer>
                    (delegate(PairwiseNucleotideSequenceAlignment alignment) { return (IPairwiseAlignmentContainer)alignment; }));
            }


        }

        /// <summary>
        /// Blast postprocessor; processes the BLAST search results;
        /// Only return the highest scoring HSP; If two or more HSP have identical scores all having this score will be returned.
        /// The number of returned HSPs will therefore occasionally exceed the number of specified best hits!!!
        /// </summary>
        public class BlastPostProcessorHighestScore : IBlastPostProcessor
        {
            private int numberOfBest = 0;
            private float highscoreAmbiguityThreshold = 10.0F;

            public BlastPostProcessorHighestScore(int numberOfBest)
            {
                this.numberOfBest = numberOfBest;
            }

            public BlastPostProcessorHighestScore(int numberOfBest, float highscoreAmbiguityThreshold)
            {
                this.numberOfBest = numberOfBest;
                this.highscoreAmbiguityThreshold = highscoreAmbiguityThreshold;
            }

            public List<IPairwiseAlignmentContainer> GetAlignments(List<PairwiseNucleotideSequenceAlignment> alignments)
            {
                if (alignments == null) return null;
                alignments.Sort(new HSPSort_ScoreDescending<PairwiseNucleotideSequenceAlignment>());

                List<IPairwiseAlignmentContainer> pwalignments = new List<IPairwiseAlignmentContainer>();

                //Determine the highscore
                float highscore = 0.0F;
                if (alignments.Count > 0) highscore = alignments[0].Score.Value;

                int n = 0;
                while (n < alignments.Count)
                {

                    if (n >= numberOfBest && (highscore - highscoreAmbiguityThreshold) > alignments[n].Score)
                    {
                        break;
                    }

                    pwalignments.Add(alignments[n]);
                    n++;
                }

                return pwalignments;
            }

        }

        /// <summary>
        /// Blast postprocessor; processes the blast search results;
        /// PairwiseNucleotideSequence alignments which are located on the same database sequence 
        /// and separated through a distance less than the specified maximum are grouped together to a CompositeNucleotideSequenceAlignment.
        /// This is to accept introns; Returns only the n-best hits; Attention the number of hits returned may exceed this specified value in case the hits with the lowest score have identical score,
        /// here all having this score are returned. If you only want the best hit always make sure that only a single hit is returned otherwise no disambigous result was obtained.
        ///</summary>
        public class BlastPostProcessor_PartialAlignmentAggregator : IBlastPostProcessor
        {
            private int? maxDistance = null;
            private int maximumOverlap = 30;
            private int best = 1;
            private float ambigiousScoreDifference = 10.0F;
            private Bio.Blast.Misc.IIntronBoundaryPolishing intronPolisher;

            /// <summary>
            /// Blast postprocessor, aggregates partial alignments to a composite alginment.
            /// Overlapping regions are resolved by keeping the highest scoring overlap.
            /// </summary>
            /// <param name="numberOfBestAlignmentsToReturn">the number of highest scoring (best) alignments to return; Attention, actually returned number may be higher depending on the ambiguity</param>
            /// <param name="substitutionMatrix">the substitution matrix which will be used to calculate the individual scores for the overlapping regions</param>
            public BlastPostProcessor_PartialAlignmentAggregator(int numberOfBestAlignmentsToReturn, SubstitutionMatrix substitutionMatrix)
            {
                intronPolisher = new IntronBoundaryPolishing_KeepBestOverlap(substitutionMatrix);
                this.best = numberOfBestAlignmentsToReturn;
            }


            /// <summary>
            /// Blast postprocessor, aggregates partial alignments to a composite alginment.
            /// Overlapping regions are resolved by keeping the highest scoring overlap.
            /// </summary>
            /// <param name="numberOfBestAlignmentsToReturn">the number of highest scoring (best) alignments to return; Attention, actually returned number may be higher depending on the ambiguity</param>
            /// <param name="substitutionMatrix">the substitution matrix which will be used to calculate the individual scores for the overlapping regions</param>
            /// <param name="maximumDistanceBetweenPartialAlignments">the maximum allowed distance between two consecutive partial alignments</param>
            public BlastPostProcessor_PartialAlignmentAggregator(int numberOfBestAlignmentsToReturn, SubstitutionMatrix substitutionMatrix, int? maximumDistanceBetweenPartialAlignments)
            {
                intronPolisher = new IntronBoundaryPolishing_KeepBestOverlap(substitutionMatrix);
                this.best = numberOfBestAlignmentsToReturn;
                this.maxDistance = maximumDistanceBetweenPartialAlignments;
            }

            /// <summary>
            /// Blast postprocessor, aggregates partial alignments to a composite alginment.
            /// Overlapping regions are resolved by keeping the highest scoring overlap.
            /// </summary>
            /// <param name="numberOfBestAlignmentsToReturn">the number of highest scoring (best) alignments to return; Attention, actually returned number may be higher depending on the ambiguity</param>
            /// <param name="substitutionMatrix">the substitution matrix which will be used to calculate the individual scores for the overlapping regions</param>
            /// <param name="maximumDistanceBetweenPartialAlignments">the maximum allowed distance between two consecutive partial alignments</param>
            /// <param name="highscoreAmgiguityDifference">minimum difference between the highest score and the subsequent ones, as to account as unambigious highscore</param>
            public BlastPostProcessor_PartialAlignmentAggregator(int numberOfBestAlignmentsToReturn, SubstitutionMatrix substitutionMatrix, int? maximumDistanceBetweenPartialAlignments, float highscoreAmgiguityDifference)
            {
                intronPolisher = new IntronBoundaryPolishing_KeepBestOverlap(substitutionMatrix);
                this.best = numberOfBestAlignmentsToReturn;
                this.maxDistance = maximumDistanceBetweenPartialAlignments;
                this.ambigiousScoreDifference = highscoreAmgiguityDifference;
            }


            /// <summary>
            /// Check if individual pairwise alignments may be parts of a higher order alignments, eg. exons of a gene.
            /// Aggregate partial alignments to a higher order composite alignment if the specified requirements are fullfilled,
            /// if not return the partial alginments unaltererd.
            /// </summary>
            /// <param name="alignments"></param>
            /// <returns></returns>
            public List<IPairwiseAlignmentContainer> GetAlignments(List<PairwiseNucleotideSequenceAlignment> alignments)
            {
                if (alignments == null) return null;
                List<IPairwiseAlignmentContainer> toRet = new List<IPairwiseAlignmentContainer>();



                //If the alignment number is 0 or 1 return the alignments unmodified;
                if (alignments.Count < 2) return alignments.ConvertAll<IPairwiseAlignmentContainer>(new Converter<PairwiseNucleotideSequenceAlignment, IPairwiseAlignmentContainer>
(delegate(PairwiseNucleotideSequenceAlignment pwal) { return (IPairwiseAlignmentContainer)pwal; }));

                List<PairwiseNucleotideSequenceAlignment> working = new List<PairwiseNucleotideSequenceAlignment>(alignments.Count);
                foreach (PairwiseNucleotideSequenceAlignment pw in alignments)
                {
                    working.Add(pw);
                }



                //Sort
                working.Sort(new HSPSort_DatabaseID_PlusPlus_DatabasePosition_QueryPosition<PairwiseNucleotideSequenceAlignment>());


                List<PairwiseNucleotideSequenceAlignment> temp = new List<PairwiseNucleotideSequenceAlignment>();


                while (working.Count > 0)//Ingenious, complicated
                {

                    temp.Add(working[0]);
                    PairwiseNucleotideSequenceAlignment last = working[0];

                    for (int i = 1; i < working.Count; i++) //Work to do
                    {
                        //DEAR GOD HAVE MERCY
                        if (last.DatabaseParent == working[i].DatabaseParent &&    //Make sure the have the same parent (either gene or chromosome)
                        last.PlusPlusStrand == working[i].PlusPlusStrand &&        //Make sure they are on the same strand, plus or minus
                        last.StartDatabase < working[i].StartDatabase &&           //make sure the have ascending start positions in the database
                        last.StartQuery < working[i].StartQuery &&                 //Make sure that the query sequence start positions are also ascending
                        last.EndQuery - maximumOverlap < working[i].StartQuery &&  //Make sure that the overlap does not exceed the maximum allowed overlap
                        last.EndDatabase < working[i].StartDatabase &&             //Make sure that the database sequence is not overlapping
                        ((maxDistance == null) || (last.EndQuery + maxDistance.Value > working[i].StartQuery)) &&           //if a maximum distance was specified make sure that it is kept in the query
                        ((maxDistance == null) || ((last.EndDatabase + maxDistance.Value) > working[i].StartDatabase)))     //if a maximum distance was specified make sure that it is kept in the database
                        {
                            //Make also sure that not both are shifted a distance exceeding the maximum distance, it would be possible to use a value other than the maximumdistance but it is more convenient this way
                            if (!((last.EndDatabase + maximumOverlap < working[i].StartDatabase) && (last.EndQuery + maximumOverlap < working[i].StartQuery)))
                            {
                                temp.Add(working[i]);
                                last = working[i];
                            }
                        }
                    }

                    if (temp.Count == 1) toRet.Add(temp[0]);
                    else
                    {
                        toRet.Add(new CompositePairwiseNucleotideSequenceAlignment(intronPolisher.GetPolishedAlignments(temp)));
                    }

                    foreach (PairwiseNucleotideSequenceAlignment pa in temp)
                    {
                        working.Remove(pa);
                    }

                    temp = new List<PairwiseNucleotideSequenceAlignment>();
                }

                toRet.Sort(new HSPSort_ScoreDescending<IPairwiseAlignmentContainer>());

                //Find the highScore;
                float highScore = 0.0F;
                if (toRet.Count > 0) highScore = toRet[0].Score.Value;

                List<IPairwiseAlignmentContainer> pwalignments = new List<IPairwiseAlignmentContainer>(best); //The number of alignments which will be returned.
                int n = 0;
                while (n < toRet.Count)
                {

                    if (n >= best && (highScore - ambigiousScoreDifference) > toRet[n].Score)
                    {
                        break;
                    }

                    pwalignments.Add(toRet[n]);
                    n++;
                }

                return pwalignments;

            }

        }


        /// <summary>
        /// Original BlastN implementation, somewhat outdated, best diagonal with the option return all would be most appropriate
        /// </summary>
        public class BlastSeedProcessor_Pairs : IBlastSeedProcessor
        {

            //SeedSorter -> use a delegate annonymous method instead


            //Seeding

            private int maximumDistanceBetweenTwoWords = 100;
            private int wordsize = 0;
            private int[] lowerWordDistanceBoundary;
            private int[] upperWordDistanceBoundary;

            public BlastSeedProcessor_Pairs(int wordSize)
            {
                this.wordsize = wordSize;
            }




            public List<BlastSeed> GetSeeds(List<BlastSeed> seeds)
            {
                //Set the word distance boundaries
                if (lowerWordDistanceBoundary == null || upperWordDistanceBoundary == null) SetWordDistanceBoundaries();
                List<BlastSeed> toReturn = new List<BlastSeed>();

                //IMPORTANT code example using an ananymous delegate for sorting
                seeds.Sort(new SortBlastSeed_DatabaseID_DatabasePosition_QueryPosition());





                //Search for neighboring seeds beeing separated through a distance of less than maximumDistance bp, and where the indel deviation is within the allowed boundaries (lower and upperBoundary)
                for (int i = 0; i < seeds.Count - 1; i++)
                {

                    int k = 1;
                    while (i + k < seeds.Count && ((seeds[i + k].DatabaseStartPosition - seeds[i].DatabaseStartPosition) < maximumDistanceBetweenTwoWords) &&
                        (seeds[i + k].DatabaseStartPosition - seeds[i].DatabaseStartPosition) >= wordsize)
                    {
                        int databaseDistance = seeds[i + k].DatabaseStartPosition - seeds[i].DatabaseStartPosition;
                        int queryDistance = seeds[i + k].QueryStartPosition - seeds[i].QueryStartPosition;


                        //Test whether the distance between two consecutive seeds is within the allowed boundaries
                        if (queryDistance >= lowerWordDistanceBoundary[databaseDistance]
                            && queryDistance <= upperWordDistanceBoundary[databaseDistance])
                        {
                            //SEED discovered!!
                            toReturn.Add(seeds[i]);


                            goto position;


                        }



                        k++;
                    }

                position: ;

                }

                return toReturn;

            }


            /// <summary>
            /// Set the boundaries for the whole range of maximum allowed distance between two consecutive words boundaries for the different words
            /// </summary>
            private void SetWordDistanceBoundaries()
            {
                this.lowerWordDistanceBoundary = new int[maximumDistanceBetweenTwoWords + 1];
                this.upperWordDistanceBoundary = new int[maximumDistanceBetweenTwoWords + 1];
                lowerWordDistanceBoundary[0] = 0;
                upperWordDistanceBoundary[0] = 0;
                for (int i = 1; i <= maximumDistanceBetweenTwoWords; i++)
                {


                    this.lowerWordDistanceBoundary[i] = (int)(i - Math.Sqrt(i));
                    this.upperWordDistanceBoundary[i] = (int)(i + Math.Sqrt(i));

                }
            }

        }




        public class BlastSeedProcessor_BestDiagonal : IBlastSeedProcessor
        {
            private int maxOffset = 5; //The allowed margin; the seed position might differ within this margin
            private int minimumDiagonal = 2;
            private int numberOfSeedsToReturn = 0;
            private bool returnAllSeeds = false;
            private float highscoreAmbigityThreshold = 1.0F;


            /// <summary>
            /// Create a BlastSeed Processor. Identifies diagonals for a given collection of BlastHashSeeds using the SSAHA algorithm. Returns a collection of BlasSeeds which may be used as basis for subsequent dynamic programing alginments algorithm
            /// Returns the 'numberofSeedsToReturn' longest identified diagonals, the number actually returned depends on the ambiguity of the diagonals.
            /// </summary>
            /// <param name="minimumDiagonalLength">minimum length of a diagonal for being identified</param>
            /// <param name="numberOfSeedsToReturn">the number of longest diagonals to return. Actually returned number may deviate from this specified value depending on diagonal ambiguity</param>
            public BlastSeedProcessor_BestDiagonal(int minimumDiagonalLength, int numberOfSeedsToReturn)
            {
                this.numberOfSeedsToReturn = numberOfSeedsToReturn;
                this.minimumDiagonal = minimumDiagonalLength;

            }

            /// <summary>
            /// Create a BlastSeed Processor. Identifies diagonals for a given collection of BlastHashSeeds using the SSAHA algorithm. Returns a collection of BlasSeeds which may be used as basis for subsequent dynamic programing alginments algorithm
            /// Returns all seeds meeting the minimum requirement in diagonal-length
            /// </summary>
            /// <param name="minimumDiagonalLength">minimum length of a diagonal for being identified</param>
            /// <param name="returnAllSeeds">specify whether all seeds should be returned. If false a default value of 0 is used</param>
            /// <param name="maxOffset">the allowed variation in base pairs of the individual seed positions within a diagonal</param>
            public BlastSeedProcessor_BestDiagonal(int minimumDiagonalLength, bool returnAllSeeds, int maxOffset)
            {
                this.returnAllSeeds = returnAllSeeds;
                this.minimumDiagonal = minimumDiagonalLength;
                this.maxOffset = maxOffset;
            }


            /// <summary>
            /// Create a BlastSeed Processor. Identifies diagonals for a given collection of BlastHashSeeds using the SSAHA algorithm. Returns a collection of BlasSeeds which may be used as basis for subsequent dynamic programing alginments algorithm
            /// Returns the 'numberofSeedsToReturn' longest identified diagonals, the number actually returned depends on the ambiguity of the diagonals.
            /// </summary>
            /// <param name="minimumDiagonalLength">minimum length of a diagonal for being identified</param>
            /// <param name="numberOfSeedsToReturn">the number of longest diagonals to return. Actually returned number may deviate from this specified value depending on diagonal ambiguity</param>
            /// <param name="maxOffset">the allowed variation in base pairs of the individual seed positions within a diagonal</param>
            public BlastSeedProcessor_BestDiagonal(int minimumDiagonalLength, int numberOfSeedsToReturn, int maxOffset)
            {
                this.numberOfSeedsToReturn = numberOfSeedsToReturn;
                this.minimumDiagonal = minimumDiagonalLength;
                this.maxOffset = maxOffset;

            }

            /// <summary>
            /// Create a BlastSeed Processor. Identifies diagonals for a given collection of BlastHashSeeds using the SSAHA algorithm. Returns a collection of BlasSeeds which may be used as basis for subsequent dynamic programing alginments algorithm
            /// Returns the 'numberofSeedsToReturn' longest identified diagonals, the number actually returned depends on the ambiguity of the diagonals.
            /// </summary>
            /// <param name="minimumDiagonalLength">minimum length of a diagonal for being identified</param>
            /// <param name="numberOfSeedsToReturn">the number of longest diagonals to return. Actually returned number may deviate from this specified value depending on diagonal ambiguity</param>
            /// <param name="maxOffset">the allowed variation in base pairs of the individual seed positions within a diagonal</param>
            /// <param name="highscoreAmbiguityThreshold">to unambigiusly identify the best seeds a minimum differenze in diagonal-length, between the best and the remaining seeds may be specified</param>
            public BlastSeedProcessor_BestDiagonal(int minimumDiagonalLength, int numberOfSeedsToReturn, int maxOffset, float highscoreAmbiguityThreshold)
            {
                this.numberOfSeedsToReturn = numberOfSeedsToReturn;
                this.minimumDiagonal = minimumDiagonalLength;
                this.maxOffset = maxOffset;
                this.highscoreAmbigityThreshold = highscoreAmbiguityThreshold;

            }

            /// <summary>
            /// Retrieve the appropriate Blast seeds. The specified number of seeds having the longest diagonal will be returned
            /// </summary>
            /// <param name="seeds">the seeds which will be investigated</param>
            /// <returns></returns>
            public List<BlastSeed> GetSeeds(List<BlastSeed> seeds)
            {
                List<BlastSeed> toReturn = new List<BlastSeed>();

                //First calculate Offsets
                foreach (BlastSeed s in seeds)
                {
                    s.Offset = s.DatabaseStartPosition - s.QueryStartPosition;
                }

                //IMPORTANT code example using an ananymous delegate for sorting
                seeds.Sort(new SortBlastSeed_DatabaseID_Offset());


                int activeDatabaseID = seeds.Count > 0 ? seeds[0].DatabaseID : 0;




                ///iterate over all Seeds find the longest diagonals and sum up diagonals which are in the same datbase sequence to a single diagonal
                for (int i = 0; i < seeds.Count - 1; i++)
                {
                    //Determine the lowest and the highest seed position
                    BlastSeed lowestDatabaseSeed = seeds[i];


                    int k = 1;
                    while (i + k < seeds.Count &&
                        seeds[i].DatabaseID == seeds[i + k].DatabaseID &&
                        seeds[i].Offset <= seeds[i + k].Offset &&
                        ((seeds[i].Offset + maxOffset) >= seeds[i + k].Offset))
                    {
                        if (seeds[i + k].DatabaseStartPosition < lowestDatabaseSeed.DatabaseStartPosition) lowestDatabaseSeed = seeds[i + k];
                        //When in preprocessormode

                        k++;
                    }

                    ///Give the first seed the weight/score of the length of the diagonal, ie the score is the number of seeds which are on the
                    ///same diagonale
                    if (k >= minimumDiagonal)
                    {
                        lowestDatabaseSeed.Score = k;

                        toReturn.Add(lowestDatabaseSeed);

                    }

                    i = i + k - 1;
                }

                toReturn.Sort(new SortBlastSeed_Score());


                ///If the number of the best which should be returned has not been set, return all BlastSeeds
                if (returnAllSeeds == true) return toReturn;
                if (numberOfSeedsToReturn == 0) return new List<BlastSeed>();

                List<BlastSeed> temp = new List<BlastSeed>();

                //A seed score difference of 1 should be allowed - necessary when using non overlapping diagonals otherwise strange results may be obtained
                float highestSeedScore = 0.0F;
                if (toReturn.Count > 0) highestSeedScore = toReturn[0].Score;

                int n = 0;
                while (n < toReturn.Count)
                {
                    if (n >= numberOfSeedsToReturn && (highestSeedScore - highscoreAmbigityThreshold) > toReturn[n].Score) //(n >= numberOfSeedsToReturn && toReturn[n - 1].Score != toReturn[n].Score)
                    {
                        break;
                    }
                    temp.Add(toReturn[n]);
                    n++;

                }
                return temp;
            }

        }

        public class BlastSeedProcessor_BestCumulativeDiagonal : IBlastSeedProcessor
        {

            /// <summary>
            /// Representation of cumulative Seeds
            /// </summary>
            private class CumulativeSeed
            {
                private List<BlastSeed> seeds;
                public CumulativeSeed()
                {
                    this.seeds = new List<BlastSeed>();
                }
                public CumulativeSeed(List<BlastSeed> seeds)
                {
                    this.seeds = seeds;
                }

                public void AddSeed(BlastSeed seed)
                {
                    seeds.Add(seed);
                }

                public float Score
                {
                    get
                    {
                        float ret = 0.0F;
                        foreach (BlastSeed s in seeds)
                        {
                            ret += s.Score;
                        }
                        return ret;
                    }
                }

                public List<BlastSeed> Seeds
                {
                    get
                    {
                        return seeds;
                    }
                }

                public BlastSeed this[int index]
                {
                    get
                    {
                        return seeds[index];
                    }
                }

                public int CountSeeds
                {
                    get
                    {
                        return seeds.Count;
                    }
                }

                /// <summary>
                /// Return the last seed
                /// </summary>
                public BlastSeed Last
                {
                    get
                    {
                        return this.seeds[seeds.Count - 1];
                    }
                }
            }

            private class CumulativeSeedSort : IComparer<CumulativeSeed>
            {
                public CumulativeSeedSort()
                {
                }

                int IComparer<CumulativeSeed>.Compare(CumulativeSeed seed1, CumulativeSeed seed2)
                {
                    //Sort in a threefold manner;
                    //First according to the databaseID;
                    //Second according to the databasePosition;
                    //and third according to the queryPosition;
                    if (seed1.Score < seed2.Score)
                    {
                        return 1;
                    }
                    else if (seed1.Score > seed2.Score)
                    {
                        return -1;
                    }
                    else
                    {
                        return 0;
                    }
                }
            }

            private int maxOffset = 5; //The allowed margin; the seed position might differ within this margin
            private int numberOfSeedsToReturn = 0;
            private int minimumCumulativeDiagonalLength; //int numberOfSeedsToReturn, int minimumDiagonalLength
            private IBlastHashtableConstructor hashConstructor;
            private float highscoreAmbigityThreshold = 1.0F;
            private int minimumDiagonalLength = 2;
            private int? maximumDistanceBetweenPartialSeeds = null;

            private BlastSeedProcessor_BestDiagonal blastSeed = null;



            /// <summary>
            /// Create a new Blast seed processor using the longest cumulative diagonal approach. The algorithm is similar to SSAHA, additionally allowing for long gaps of the specified maximum length;
            /// </summary>
            /// <param name="minimumCumulativeDiagonalLength">minimum number of words for a cumulative diagonal</param>
            /// <param name="numberOfSeedsToReturn">the number of cumulative diagonals which should be returned.
            /// Attention cumulative diagonals are split into diagonals when returned, the actually returned number of seeds might thus be much larger than the specified value.
            /// This diagonals may be used as individual seeds for dynamic programing algorithm.
            /// </param>
            /// <param name="hashConstructor">the class used for constructing the Blast hash table (word table). This is necessary since this algorithm requires information about the wordlength and whether the words are overlapping</param>
            /// <param name="maximumDistanceBetweenPartialDiagonals">the maximum distance between two isolated diagonals as to be considered for aggregation into a cumulative diagonal; if null an indefinitely large distance will be allowed</param>
            public BlastSeedProcessor_BestCumulativeDiagonal(int minimumCumulativeDiagonalLength, int numberOfSeedsToReturn, IBlastHashtableConstructor hashConstructor, int? maximumDistanceBetweenPartialDiagonals)
            {
                this.maximumDistanceBetweenPartialSeeds = maximumDistanceBetweenPartialDiagonals;
                this.numberOfSeedsToReturn = numberOfSeedsToReturn;
                this.minimumCumulativeDiagonalLength = minimumCumulativeDiagonalLength;
                this.hashConstructor = hashConstructor;

                //Initialize the seed pre-processing phase, the null specifies that the pre-processor will return all seeds meeting the minimum required diagonal length


            }


            /// <summary>
            /// Create a new Blast seed processor using the longest cumulative diagonal approach. The algorithm is similar to SSAHA, additionally allowing for long gaps of the specified maximum length;
            /// </summary>
            /// <param name="minimumCumulativeDiagonalLength">minimum number of words for a cumulative diagonal</param>
            /// <param name="numberOfSeedsToReturn">the number of cumulative diagonals which should be returned.
            /// Attention cumulative diagonals are split into diagonals when returned, the actually returned number of seeds might thus be much larger than the specified value.
            /// This diagonals may be used as individual seeds for dynamic programing algorithm.
            /// </param>
            /// <param name="hashConstructor">the class used for constructing the Blast hash table (word table). This is necessary since this algorithm requires information about the wordlength and whether the words are overlapping</param>
            /// <param name="maximumDistanceBetweenPartialDiagonals">the maximum distance between two isolated diagonals as to be considered for aggregation into a cumulative diagonal; if null an indefinitely large distance will be allowed</param>
            /// <param name="maxOffset">maxium deviation of the individual words from a perfect diagonal</param>
            public BlastSeedProcessor_BestCumulativeDiagonal(int minimumCumulativeDiagonalLength, int numberOfSeedsToReturn, IBlastHashtableConstructor hashConstructor, int? maximumDistanceBetweenPartialDiagonals, int maxOffset)
            {

                this.maximumDistanceBetweenPartialSeeds = maximumDistanceBetweenPartialDiagonals;
                this.numberOfSeedsToReturn = numberOfSeedsToReturn;
                this.minimumCumulativeDiagonalLength = minimumCumulativeDiagonalLength;
                this.hashConstructor = hashConstructor;
                this.maxOffset = maxOffset;
            }


            /// <summary>
            /// Create a new Blast seed processor using the longest cumulative diagonal approach. The algorithm is similar to SSAHA, additionally allowing for long gaps of the specified maximum length;
            /// </summary>
            /// <param name="minimumCumulativeDiagonalLength">minimum number of words for a cumulative diagonal</param>
            /// <param name="numberOfSeedsToReturn">the number of cumulative diagonals which should be returned.
            /// Attention cumulative diagonals are split into diagonals when returned, the actually returned number of seeds might thus be much larger than the specified value.
            /// This diagonals may be used as individual seeds for dynamic programing algorithm.
            /// </param>
            /// <param name="hashConstructor">the class used for constructing the Blast hash table (word table). This is necessary since this algorithm requires information about the wordlength and whether the words are overlapping</param>
            /// <param name="maximumDistanceBetweenPartialDiagonals">the maximum distance between two isolated diagonals as to be considered for aggregation into a cumulative diagonal; if null an indefinitely large distance will be allowed</param>
            /// <param name="maxOffset">maxium deviation of the individual words from a perfect diagonal</param>
            /// <param name="highScoreAmbiguityThreshold">the minimum difference between the highest score and the remaining scores as to account as unambigious</param>
            public BlastSeedProcessor_BestCumulativeDiagonal(int minimumCumulativeDiagonalLength, int numberOfSeedsToReturn, IBlastHashtableConstructor hashConstructor, int? maximumDistanceBetweenPartialDiagonals, int maxOffset, int highScoreAmbiguityThreshold)
            {

                this.maximumDistanceBetweenPartialSeeds = maximumDistanceBetweenPartialDiagonals;
                this.numberOfSeedsToReturn = numberOfSeedsToReturn;
                this.minimumCumulativeDiagonalLength = minimumCumulativeDiagonalLength;
                this.hashConstructor = hashConstructor;
                this.maxOffset = maxOffset;
                this.highscoreAmbigityThreshold = highScoreAmbiguityThreshold;
            }



            /// <summary>
            /// Return the specified number of Blast seeds
            /// </summary>
            /// <param name="seeds"></param>
            /// <returns></returns>
            public List<BlastSeed> GetSeeds(List<BlastSeed> seeds)
            {
                //Initialize the preprozessor if not yet initialised
                if (this.blastSeed == null) blastSeed = new BlastSeedProcessor_BestDiagonal(this.minimumDiagonalLength, true, this.maxOffset);

                List<CumulativeSeed> cumSeed = new List<CumulativeSeed>();
                List<BlastSeed> working = this.blastSeed.GetSeeds(seeds);

                ///Initiate the cumulative seed ends
                foreach (BlastSeed bs in working)
                {
                    bs.DatabaseEndPosition = bs.DatabaseStartPosition + ((int)bs.Score - 1) * hashConstructor.Stepsize + hashConstructor.Wordsize - 1;
                    bs.QueryEndPosition = bs.QueryStartPosition + ((int)bs.Score - 1) * hashConstructor.Stepsize + hashConstructor.Wordsize - 1;
                }
                working.Sort(new SortBlastSeed_DatabaseID_DatabasePosition_QueryPosition());

                if (working.Count == 0) return working;


#if DEBUG           //Test whether if seed consists of a diagonal
                foreach (BlastSeed bs in working)
                {
                    if (bs.Score < 2) throw new InvalidOperationException("Each seed should consist of at least of one diagonal");
                }
#endif


                List<BlastSeed> temp = new List<BlastSeed>();
                while (working.Count > 0)
                {
                    temp.Add(working[0]);
                    BlastSeed last = working[0];
                    for (int i = 1; i < working.Count; i++)
                    {
                        if ((last.DatabaseID == working[i].DatabaseID) &&
                            (last.DatabaseEndPosition < working[i].DatabaseStartPosition) &&
                            (last.QueryEndPosition < working[i].QueryStartPosition)
                            && ((this.maximumDistanceBetweenPartialSeeds == null) ||
                            (last.DatabaseEndPosition + this.maximumDistanceBetweenPartialSeeds.Value > working[i].DatabaseStartPosition)))
                        {

                            last = working[i];
                            temp.Add(last);
                        }

                    }
                    cumSeed.Add(new CumulativeSeed(temp));

                    foreach (BlastSeed s in temp)
                    {
                        working.Remove(s);
                    }
                    temp = new List<BlastSeed>();


                }


                /*
                for (int i = 1; i < working.Count; i++)
                {
                    if((active.Last.DatabaseID == working[i].DatabaseID) &&
                        (active.Last.DatabaseEndPosition < working[i].DatabaseStartPosition) &&
                        (active.Last.QueryEndPosition < working[i].QueryStartPosition) &&
                    ((this.maximumDistanceBetweenPartialSeeds==null)||
                    (active.Last.DatabaseEndPosition+this.maximumDistanceBetweenPartialSeeds>working[i].DatabaseStartPosition)))
                    {
                        active.AddSeed(working[i]);
                    }
                    else
                    {
                        if (active.Score >= this.minimumCumulativeDiagonalLength) cumSeed.Add(active); //If the seed meets the minimum requirements add it to the list
                        active = new CumulativeSeed();
                        active.AddSeed(working[i]);
                    }

                }

                 */



                List<CumulativeSeed> tempCumulative = new List<CumulativeSeed>();
                foreach (CumulativeSeed cs in cumSeed)
                {
                    if (cs.Score >= minimumCumulativeDiagonalLength) tempCumulative.Add(cs);
                }

                ///Sort the cumulative Seeds
                tempCumulative.Sort(new CumulativeSeedSort());



                ///If the number of the best which should be returned has not been set, return all BlastSeeds
                if (numberOfSeedsToReturn == 0) throw new InvalidOperationException("the number of the best seeds which should be returned has to be set to a value");

                temp = new List<BlastSeed>();
                float highestSeedScore = 0.0F;
                if (tempCumulative.Count > 0) highestSeedScore = tempCumulative[0].Score;

                ///Return the specified number of best cumulative seeds.
                ///Directly before returning them to the BLAST engine they will be split into isolated seeds.
                int n = 0;
                while (n < tempCumulative.Count)
                {

                    if (n >= numberOfSeedsToReturn && (highestSeedScore - highscoreAmbigityThreshold) > tempCumulative[n].Score)
                    {
                        break;
                    }
                    temp.AddRange(tempCumulative[n].Seeds);
                    n++; //The best seeds have to be counted independent of the variable temp.Count since this will be higher

                }

                return temp;


            }

        }


        public interface IAnchoredDynamicProgramming : IDynamicProgramming
        {
            IAnchoredDynamicProgramming NewDynamicProgramming(ISequenceContainer database, ISequenceContainer query, SubstitutionMatrix substitutionMatrix, int anchorPositionDatabase, int anchorPositionQuery);
            IAnchoredDynamicProgramming NewDynamicProgramming(ISequenceContainer database, ISequenceContainer query, SubstitutionMatrix substitutionMatrix, int anchorPositionDatabase, int anchorPositionQuery, bool plusplusOrientation);
        }


 

}