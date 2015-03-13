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
using System.Text.RegularExpressions;

namespace Bio.Seq.Statistics
{


 
        /// <summary>
        /// Provides basic functionallity for grouping together similar nucleotide sequences,
        /// calculating mean values, standard deviations thereof, the density [counts/Mbp] are calculated
        /// New nucleotide sequences are added with the .Add() method, the properties might be use to retrieve the statistics
        /// </summary>
        public abstract class PositionableStatistic_Base<TPositionable>
            where TPositionable : IPositionable
        {
            private List<TPositionable> nucleotides = new List<TPositionable>();
            private INucSeqInfo nucSeqInfo = null;
            private long? totalLength = null;

            private double? averageLength = null;
            private double? averageLengthStandardDeviationVariable = null;
            private double? averageLengthStandardDeviationSample = null;


            private double? densityAllNucleotides = null;
            private double? densityWithoutN = null;
            private double? densityATUCG = null;



            public PositionableStatistic_Base(List<TPositionable> toInvestigate)
            {
                this.nucleotides = toInvestigate;
            }
            public PositionableStatistic_Base(List<TPositionable> toInvestigate, INucSeqInfo nucSeqInfo)
                : this(toInvestigate)
            {
                this.nucSeqInfo = nucSeqInfo;
            }


            public PositionableStatistic_Base()
            {
            }








            #region Properties

            public INucSeqInfo NucSeqInfo
            {
                get
                {
                    return this.nucSeqInfo;
                }
                protected set
                {
                    this.nucSeqInfo = value;
                }
            }

            /// <summary>
            /// Calcualates the total length of all nucleotides in the collection taken together
            /// Lazy initiated and stored as variable- not calculated anew with each request.
            /// </summary>
            public long TotalLengthOfAllNucleotides
            {
                get
                {
                    if (totalLength == null)
                    {
                        this.totalLength = 0;
                        foreach (TPositionable ns in nucleotides)
                        {
                            totalLength += ns.Length;
                        }
                    }
                    return totalLength.Value;

                }
            }


            /// <summary>
            /// Returns the actual collection of nucleotide sequences contained in this class as read only!
            /// </summary>
            public List<TPositionable> NucleotideSequences
            {
                get
                {
                    return nucleotides;
                }

            }




            /// <summary>
            /// Indexer: provides access to the elements of the statistics collection
            /// </summary>
            /// <param name="index"></param>
            /// <returns></returns>
            public TPositionable this[int index]
            {

                get
                {
                    return nucleotides[index];
                }

            }

            /// <summary>
            /// The average length of all nucleotide sequences in this collection
            /// Lazy initated and stored - not calculated anew with each request
            /// </summary>
            public double? AverageLength
            {
                get
                {
                    if (nucleotides.Count == 0) return null;
                    if (this.averageLength == null)
                    {
                        this.averageLength = ((double)this.TotalLengthOfAllNucleotides) / (double)nucleotides.Count;
                    }

                    return this.averageLength;
                }
            }

            /// <summary>
            /// Returns the standard deviation of the Average length; 1/N
            /// Assuming that the collection represents the whole population
            /// <see cref="http://en.wikipedia.org/wiki/Standard_deviation"/>
            /// </summary>
            public double? StdDev_Length_RandomVariable
            {
                get
                {

                    if (this.AverageLength == null) return null;
                    if (this.averageLengthStandardDeviationVariable == null)
                    {
                        double devSum = 0.0;

                        foreach (TPositionable nucSeq in nucleotides)
                        {
                            devSum += Math.Pow((((double)nucSeq.Length) - (double)averageLength), 2.0);
                        }
                        this.averageLengthStandardDeviationVariable = Math.Sqrt((devSum / ((double)this.Count)));
                    }

                    return this.averageLengthStandardDeviationVariable;
                }
            }


            /// <summary>
            /// Returns the standard deviation of the Average length; 1/(N-1)
            /// Assuming that the collection represents only a sample of a larger population
            /// <see cref="http://en.wikipedia.org/wiki/Standard_deviation"/>
            /// </summary>
            public double? StdDev_Length_SampleOfPopulation
            {
                get
                {

                    if (this.AverageLength == null) return null;
                    if (this.averageLengthStandardDeviationSample == null)
                    {
                        double devSum = 0.0;

                        foreach (TPositionable nucSeq in nucleotides)
                        {
                            devSum += Math.Pow((((double)nucSeq.Length) - (double)averageLength), 2.0);
                        }
                        this.averageLengthStandardDeviationSample = Math.Sqrt((devSum / (((double)this.Count) - 1.0)));
                    }

                    return this.averageLengthStandardDeviationSample;
                }
            }

            /// <summary>
            /// Returns the number of nucleotide sequence elements in this class
            /// </summary>
            public int Count
            {
                get
                {
                    return nucleotides.Count;
                }
            }


            /// <summary>
            /// Frequency of a INucleotideSequence element per Mbp
            /// not considering the N nucleotides
            /// </summary>
            public double? DensityAllNucleotides
            {
                get
                {
                    if (this.nucSeqInfo == null || this.nucSeqInfo.Length == 0) return null;
                    if (densityAllNucleotides == null)
                    {
                        densityAllNucleotides = (1000000.0 / ((double)this.nucSeqInfo.Length / (double)this.nucleotides.Count));
                    }
                    return densityAllNucleotides;
                }
            }

            /// <summary>
            /// Frequency of a INucleotideSequence element per Mbp
            /// not considering the N nucleotides
            /// </summary>
            public double? DensityWithoutN
            {
                get
                {
                    if (this.nucSeqInfo == null || this.nucSeqInfo.Count_AllChar == 0) return null;
                    if (this.densityWithoutN == null)
                    {
                        densityWithoutN = (1000000.0 / ((double)(this.nucSeqInfo.Count_AllChar - this.nucSeqInfo.Count_N) / (double)this.nucleotides.Count));
                    }
                    return densityWithoutN;
                }
            }

            /// <summary>
            /// Frequency of a INucleotideSequence element per Mbp
            /// only considering ATUCG nucleotides
            /// </summary>
            public double? DensityATUCG
            {
                get
                {
                    if (this.nucSeqInfo == null || this.nucSeqInfo.Count_ATUCG == 0) return null;
                    if (this.densityATUCG == null)
                    {
                        densityATUCG = (1000000.0 / ((double)(this.nucSeqInfo.Count_ATUCG) / (double)this.nucleotides.Count));
                    }
                    return densityATUCG;
                }
            }


            /// <summary>
            /// An element of this NucleotideStatistic collection is on the average occuring all .... base pairs
            /// </summary>
            public double? NucleotidesPerCount
            {
                get
                {
                    if (this.nucSeqInfo == null || this.nucSeqInfo.Length == 0) return null;
                    return ((double)this.nucSeqInfo.Length / (double)this.nucleotides.Count);

                }
            }

            /// <summary>
            /// An element of this NucleotideStatistic collection is on the average occuring all ....nucleotide sequences. (e.g.: single fastas)
            /// </summary>
            public double? FilesPerCount
            {
                get
                {
                    if (nucleotides.Count == 0 || nucSeqInfo == null) return null;
                    return ((double)nucSeqInfo.SequenceCount / (double)nucleotides.Count);

                }
            }



            #endregion;
        }



        /// <summary>
        /// Simple implementation of the Nucleotide sequence statistic class
        /// </summary>
        /// <typeparam name="TNucleotideSequence"></typeparam>
        public class NucleotideSequenceStatistic<TNucleotideSequence> : PositionableStatistic_Base<TNucleotideSequence> where TNucleotideSequence : INucleotideSequence
        {
            public double? averagegcContent = null;
            public NucleotideSequenceStatistic(List<TNucleotideSequence> toInvestigate)
                : base(toInvestigate)
            {
            }
            public NucleotideSequenceStatistic(List<TNucleotideSequence> toInvestigate, INucSeqInfo nucSeqInfo)
                : base(toInvestigate, nucSeqInfo)
            {
            }
            public NucleotideSequenceStatistic()
            {
            }

            public double AverageGC_Content
            {
                get
                {
                    if (averagegcContent == null)
                    {
                        double gctotal = 0;
                        foreach (TNucleotideSequence ts in this.NucleotideSequences)
                        {
                            gctotal += ts.Count_CG;
                            this.averagegcContent = gctotal / this.TotalLengthOfAllNucleotides;
                        }
                    }
                    return averagegcContent.Value;
                }
            }
        }


        /// <summary>
        /// Provides statistics for the distribution of IPositionable elements within a chromosome or more generally any large sequence.
        /// Holds a collection of ChromosomeBin elements which provide the actual information such as the average length, frequency etc of the IPositionables in the sequence
        /// </summary>
        /// <typeparam name="TPositionable"></typeparam>
        public class ChromosomeDistributionStatistic<TPositionable> : PositionableStatistic_Base<TPositionable>, IPositionable where TPositionable : IPositionable
        {

            private INucleotideSequence nucleotideSequence;
            private List<ChromosomeBin<TPositionable>> bins = null;
            private int? startPosition = null;
            private int? endPosition = null;
            //private string parentName=null;
            private bool? isRoot = null;


            public ChromosomeDistributionStatistic(List<ChromosomeBin<TPositionable>> bins, List<TPositionable> allPositionables)
                : base(allPositionables)
            {
                this.bins = bins;
                bins.Sort(new Sort.SortPositionables_ParentName_StartPos_Ascending<ChromosomeBin<TPositionable>>());

                allPositionables.Sort(new Sort.SortPositionables_ParentName_StartPos_Ascending<TPositionable>());

                if (allPositionables.Count > 0)
                {
                    this.NucSeqInfo = new NucSeqInfo(allPositionables.Count > 0 ? allPositionables[0].ParentName : "");
                    this.NucSeqInfo.Length = allPositionables[allPositionables.Count - 1].End.Value;
                }
                InitializeAllBins();
            }




            private void InitializeAllBins()
            {
                foreach (ChromosomeBin<TPositionable> cbin in bins)
                {
                    cbin.ChromosomeDistributionStatistic_Reference = this;
                }
            }





            /// <summary>
            /// The name of the chromosome sequence,
            /// that is the fasta identifier
            /// </summary>
            public string ChromosomeName
            {
                get
                {
                    if (this.NucSeqInfo == null) return null;


                    return this.NucSeqInfo.Name;

                }
                set
                {
                    if (!(this.NucSeqInfo == null))
                    {
                        this.NucSeqInfo.Name = value;
                    }


                }
            }


            /// <summary>
            /// The whole sequence of the parent sequence,
            /// that is the whole sequence of the chromosome
            /// </summary>
            public INucleotideSequence NucleotideSequence
            {
                get
                {
                    return this.nucleotideSequence;
                }
                set
                {

                    if (value.Sequence.Length > 0)
                    {
                        this.nucleotideSequence = value;
                    }
                }
            }


            public new ChromosomeBin<TPositionable> this[int index]
            {
                get
                {
                    return this.bins[index];
                }
            }

            /// <summary>
            /// Length of the chromosome or sequence
            /// </summary>
            public long Length
            {
                get
                {
                    if (this.NucSeqInfo == null) return 0;
                    return NucSeqInfo.Length;
                }

            }

            /// <summary>
            /// 
            /// </summary>
            public int? Start
            {
                get
                {
                    return this.startPosition;
                }
                set
                {
                    this.startPosition = value;
                }
            }

            public int? End
            {
                get
                {
                    return this.endPosition;
                }
                set
                {
                    this.endPosition = value;
                }
            }

            /// <summary>
            /// The parent name;
            /// use only when this is not a whole chromosome but a subsequence of a larger sequence
            /// </summary>
            public string ParentName
            {
                get
                {
                    return this.NucSeqInfo.ParentName;
                }
                set
                {
                    this.NucSeqInfo.ParentName = value;
                }
            }
            public bool? IsRoot
            {
                get
                {
                    return this.isRoot;
                }
                set
                {
                    this.isRoot = value;
                }
            }

            /// <summary>
            /// The number of ChromosomalBins upon which this statistic is based
            /// </summary>
            public new int Count
            {
                get
                {
                    return bins.Count;
                }
            }




      



            public string FeatureName
            {
                get { return "Chromosome Distribution"; }
            }
        }






        /// <summary>
        /// Provides basic statistics for a given set of nucleotide sequences being part of a bin.
        /// A bin is a physically part of a chromosome or nucleotide sequence. 
        /// For example with a bin size of 1Mbp the human chromosome 1 can be divided into at least 220 bins of the size 1Mbp
        /// Since the step size of the bin is not automatically equal to the bin size the number of bins can be much larger, that is bins might be overlapping
        /// -> slidding window approach
        /// </summary>
        /// <typeparam name="TNucleotidSequence">implements the interface INucleotideSequence</typeparam>
        public class ChromosomeBin<TPositionable> : PositionableStatistic_Base<TPositionable>, IPositionable
            where TPositionable : IPositionable
        {
            private int? end = null;
            private int? start = null;
            private ChromosomeDistributionStatistic<TPositionable> referenceChromosomeStatistic;

            public ChromosomeBin(List<TPositionable> toInvestigate, string parentName, int start, int length)
                : base(toInvestigate)
            {



                this.NucSeqInfo = new NucSeqInfo("");
                this.NucSeqInfo.ParentName = parentName;
                this.NucSeqInfo.Length = length;

                this.end = start + length - 1;
                this.start = start;

            }

            public ChromosomeBin(List<TPositionable> toInvestigate, ChromosomeDistributionStatistic<TPositionable> referenceChromosomeStatistic, int start, int length)
                : this(toInvestigate, referenceChromosomeStatistic.ParentName, start, length)
            {
                this.referenceChromosomeStatistic = referenceChromosomeStatistic;
            }




            public ChromosomeDistributionStatistic<TPositionable> ChromosomeDistributionStatistic_Reference
            {
                get
                {
                    return this.referenceChromosomeStatistic;
                }
                set
                {
                    this.referenceChromosomeStatistic = value;
                }
            }






            public int? Start
            {
                get
                {
                    return this.start;
                }
                set
                {
                    this.start = value;
                }
            }

            public int? End
            {
                get
                {
                    return this.end;
                }
                set
                {
                    this.end = value;
                }
            }

            /// <summary>
            /// The middle position of the bin
            /// </summary>
            public double? PositionMiddle
            {
                get
                {
                    if (this.start == null || this.end == null) return null;
                    return ((double)this.end + ((double)this.start)) / 2.0;
                }

            }

            public bool? IsRoot
            {
                get
                {
                    return false;
                }
                set
                {
                }
            }

            public string ParentName
            {
                get
                {
                    return this.NucSeqInfo.ParentName;
                }
                set
                {
                    this.NucSeqInfo.ParentName = value;
                }
            }



            /// <summary>
            /// the length of the chromosomal bin
            /// </summary>
            public long Length
            {

                get
                {

                    return this.NucSeqInfo.Length;
                }

            }





  


            public string FeatureName
            {
                get { return "Chromosome Bin"; }
            }

  
        }





        /// <summary>
        /// Creates physical bin statistics
        /// </summary>
        public class ChromosomeStatisticProvider<TPositionable>
            where TPositionable : IPositionable
        {

            private List<TPositionable> toInvestigate;


            public ChromosomeStatisticProvider(List<TPositionable> toInvestigate)
            {

                this.toInvestigate = toInvestigate;


            }




            /// <summary>
            /// 
            /// </summary>
            /// <returns>Bin statistic, List.Count=0 if no more BinStatistics are available</returns>
            public ChromosomeDistributionStatistic<TPositionable> GetBinStatistic(string sequenceToInvestigate, int binSize, int stepSize)
            {
                List<TPositionable> chromosomalSequences = new List<TPositionable>();

                foreach (TPositionable seq in toInvestigate)
                {
                    if (seq.ParentName == sequenceToInvestigate) chromosomalSequences.Add(seq);
                }

                   chromosomalSequences.Sort(new Sort.SortPositionables_ParentName_StartPos_Ascending<TPositionable>());

                List<ChromosomeBin<TPositionable>> bins = new List<ChromosomeBin<TPositionable>>();





                int binStart = 1;
                int binEnd = 1 + binSize;
                int nextBinStart = 1 + stepSize;
                int nextBinStartIndex = 0;

                List<TPositionable> binSequences = new List<TPositionable>();

                for (int i = 0; i < chromosomalSequences.Count; i++)
                {

                    if (chromosomalSequences[i].Start >= binStart && chromosomalSequences[i].Start < binEnd)
                    {
                        binSequences.Add(chromosomalSequences[i]);
                    }
                    else if (chromosomalSequences[i].Start >= binEnd)
                    {


                        bins.Add(new ChromosomeBin<TPositionable>(binSequences, sequenceToInvestigate, binStart, binSize));

                        binSequences = new List<TPositionable>();

                        nextBinStart += stepSize;
                        binStart += stepSize;
                        binEnd += stepSize;

                        i = nextBinStartIndex;

                        continue;

                    }

                    if (chromosomalSequences[i].Start < nextBinStart) nextBinStartIndex = i;

                }
                //Add the last element, here the bin size is not the normal size

                if (binSequences.Count > 0) bins.Add(new ChromosomeBin<TPositionable>(binSequences, sequenceToInvestigate, binStart, (int)(binSequences[binSequences.Count - 1].End - binStart + 1)));
                binSequences = null;



                return new ChromosomeDistributionStatistic<TPositionable>(bins, chromosomalSequences);

            }



        }



    

}