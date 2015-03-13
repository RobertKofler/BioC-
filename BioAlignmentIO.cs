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
using System.IO;


namespace Bio.Alignment.IO
{
    using Bio.Seq;
    using Bio.Alignment;
    using Bio.Seq.IO;


    /// <summary>
    /// Writes the alignment of a nucleotide sequence in the tab delimited format
    /// </summary>
    public class PairwiseNucleotideSequenceAlignmentWriter
    {
        private StreamWriter sw;
        private bool wroteIni = false;
        public PairwiseNucleotideSequenceAlignmentWriter(string path) : this(new StreamWriter(path, false, Encoding.ASCII)) { }
        public PairwiseNucleotideSequenceAlignmentWriter(StreamWriter sw)
        {
            this.sw = sw;
        }


        public void WriteAlignments(List<IPairwiseNucleotideSequenceAlignmentContainer> alignments)
        {

#if DEBUG
            System.Diagnostics.Stopwatch stopw = new System.Diagnostics.Stopwatch();
            stopw.Start();
#endif
            if (!wroteIni)
            {
                sw.WriteLine("#>> ...Start of alignment");
                sw.WriteLine("#dbName: ...database sequence name");
                sw.WriteLine("#queName: ...query sequence name");
                sw.WriteLine("#subAlgn: ...number of sub-alignments (>0)");
                sw.WriteLine("#strand: ...strand, either plus or minus");
                sw.WriteLine("#dbLength: ...total length of the database sequence");
                sw.WriteLine("#queLength: ...total length of the whole query sequence");
                sw.WriteLine("#summary:\tdatabase_start\tdatabase_end\tquery_start\tquery_end\tscore\tsimilarity\tgaps");
                sw.WriteLine("#>\tdatabase_start\tdatabase_end\tquery_start\tquery_end\tscore\tsimilarity\tgaps ...start of sub-alignment");
                sw.WriteLine("#NNNNNNNNNNNNN ...sub-alignment database sequence");
                sw.WriteLine("#||||||||||||| ...sub-alignment similarity sequence");
                sw.WriteLine("#NNNNNNNNNNNNN ...sub-alignment query sequence");
                sw.WriteLine("#< ... end of subalignment");
                sw.WriteLine("#<< ...end of alignment");
                sw.WriteLine();
                sw.WriteLine();
                wroteIni = true;
            }


            foreach (IPairwiseNucleotideSequenceAlignmentContainer al in alignments)
            {
                WriteAlignment(al);
            }

#if DEBUG
            stopw.Stop();
            System.Diagnostics.Debug.WriteLine(String.Format("{0} FastaWriter; {1} files saved in {2} ms; First file contained {3} nucleotides; Sequence used: {4}", DateTime.Now.TimeOfDay, alignments.Count, stopw.ElapsedMilliseconds,
                alignments.Count <= 0 ? 0 : alignments[0].Length, alignments.Count <= 0 ? "no sequence" : alignments[0].Alignment.GetType().ToString()));


#endif
        }

        public void WriteAlignment(IPairwiseNucleotideSequenceAlignmentContainer alignment)
        {
            sw.WriteLine(">>");
            sw.WriteLine(String.Format("dbName: {0}", alignment.DatabaseParent));
            sw.WriteLine(String.Format("queName: {0}", alignment.QueryParent));
            sw.WriteLine(String.Format("subAlgn: {0}", alignment.Count_SubAlignments));
            sw.WriteLine(String.Format("strand: {0}", alignment.PlusPlusStrand == true ? "+/+" : "+/-"));
            if (alignment.LengthDatabaseParent != null) sw.WriteLine(String.Format("dbLength: {0}", alignment.LengthDatabaseParent));
            if (alignment.LengthQueryParent != null) sw.WriteLine(String.Format("queLength: {0}", alignment.LengthQueryParent));
            sw.WriteLine(String.Format("summary:\t{0}\t{1}\t{2}\t{3}\t{4:f}\t{5:f}\t{6}", alignment.StartDatabase, alignment.EndDatabase, alignment.StartQuery, alignment.EndQuery, alignment.Score, alignment.SimilarityWithoutGaps, alignment.Gaps));
            for (int i = 0; i < alignment.Count_SubAlignments; i++)
            {
                WriteSubAlignment(alignment[i]);
            }
            sw.WriteLine("<<");
            sw.WriteLine();
            sw.WriteLine();
        }

        private void WriteSubAlignment(IPairwiseAlignmentContainer pw)
        {
            int length = 70;
            ISequenceContainer dat = pw.Alignment.DatabaseSequence;
            ISequenceContainer que = pw.Alignment.QuerySequence;
            ISequenceContainer sim = pw.Alignment.SimilaritySequence;

            int count = 0;
            sw.WriteLine(String.Format(">\t{0}\t{1}\t{2}\t{3}\t{4:f}\t{5:f}\t{6}", pw.StartDatabase, pw.EndDatabase, pw.StartQuery, pw.EndQuery, pw.Score, pw.SimilarityWithoutGaps, pw.Gaps));
            while (dat.Length >= count + length)
            {

                sw.WriteLine(dat.ToString(count, length));
                sw.WriteLine(sim.ToString(count, length));
                sw.WriteLine(que.ToString(count, length));
                sw.WriteLine();

                count += length;

            }

            if (dat.Length - count > 0)
            {
                sw.WriteLine(dat.ToString(count, dat.Length - count));
                sw.WriteLine(sim.ToString(count, sim.Length - count));
                sw.WriteLine(que.ToString(count, que.Length - count));
            }
            sw.WriteLine("<");
        }


    }


    /// <summary>
    /// Reader for NucleotideSequencePairwiseAlignments
    /// Reades composite and single nucleotide sequence alignments
    /// </summary>
    public class PairwiseNucleotideSequenceAlignmentReader
    {

        /// <summary>
        /// Builder for PairwiseNucleotideSequenceAlignments
        /// </summary>
        private class SubAlignmentBuilder
        {
            private ISequenceContainer datab;
            private ISequenceContainer query;
            private string databaseParent;
            private string queryParent;
            private bool plusplus;
            private int databaseStart;
            private int databaseEnd;
            private int queryStart;
            private int queryEnd;
            private int? lengthDatabaseParent;
            private int? lengthQueryParent;
            private float score;



            public SubAlignmentBuilder(string databaseParent, string queryParent, bool plusplus, int? lengthDatabaseParent, int? lengthQueryParent, string initial)
            {//>	909	950	1	42	126,00	100,00	0

                this.databaseParent = databaseParent;
                this.queryParent = queryParent;
                this.plusplus = plusplus;
                this.lengthDatabaseParent = lengthDatabaseParent;
                this.lengthQueryParent = lengthQueryParent;
                string[] frags = initial.Split(new char[] { '\t' });
                this.databaseStart = Convert.ToInt32(frags[1]);
                this.databaseEnd = Convert.ToInt32(frags[2]);
                this.queryStart = Convert.ToInt32(frags[3]);
                this.queryEnd = Convert.ToInt32(frags[4]);
                this.score = (float)Convert.ToDouble(frags[5]);
                this.datab = SequenceFactory.GetDefaultSequence();
                this.query = SequenceFactory.GetDefaultSequence();

            }
            public SubAlignmentBuilder()
            {
            }

            public void AddDatabasesequence(string databaseSequence)
            {
                this.datab.Append(databaseSequence);
            }
            public void AddQuerysequence(string querySequence)
            {
                this.query.Append(querySequence);
            }
            public PairwiseNucleotideSequenceAlignment GetAlignment()
            {
                PairwiseNucleotideSequenceAlignment pw = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(datab, query), databaseParent, queryParent, databaseStart, queryStart, databaseEnd, queryEnd);
                pw.PlusPlusStrand = plusplus;
                pw.Score = score;
                pw.LengthDatabaseParent = lengthDatabaseParent;
                pw.LengthQueryParent = lengthQueryParent;
                return pw;
            }
        }

        /// <summary>
        /// Builder for Composite or PairwiseAlignments depending on the number of sub-alignments
        /// 1..PairwiseNucleotideSequence alignment
        /// >1..CompositeNucleotideSequence alignment
        /// </summary>
        private class AlignmentContainerBuilder
        {
            private List<PairwiseNucleotideSequenceAlignment> pwals;
            public AlignmentContainerBuilder()
            {
                pwals = new List<PairwiseNucleotideSequenceAlignment>();
            }

            /// <summary>
            /// Add another Subalignment
            /// </summary>
            /// <param name="pw"></param>
            public void AddSubalignment(PairwiseNucleotideSequenceAlignment pw)
            {
                pwals.Add(pw);
            }

            /// <summary>
            /// Retrieve the appropriate PairwiseAlignment
            /// either Composite or PairwiseNucleotideSequenceAlignment
            /// </summary>
            /// <returns></returns>
            public IPairwiseAlignmentContainer GetPairwiseAlignmentContainer()
            {
                if (pwals.Count == 1) return pwals[0];
                else if (pwals.Count > 1) return new CompositePairwiseNucleotideSequenceAlignment(pwals);
                else throw new InvalidOperationException("Invalid number of subalignments");
            }
        }

        /// <summary>
        /// Working variables
        /// </summary>
        private enum ReadMode { Nonsense, GenerallInfo, SubAlignment };
        private StreamReader sr;

        private readonly int bufferSize;

        /// <summary>
        /// Create new NucleotideSequenceAlignmentReader
        /// </summary>
        /// <param name="sr">a StreamReader instance</param>
        /// <param name="bufferSize">the bufferSize</param>
        public PairwiseNucleotideSequenceAlignmentReader(StreamReader sr, int bufferSize)
        {
            this.sr = sr;
            this.bufferSize = bufferSize;

        }

        /// <summary>
        /// Create new NucleotideSequenceAlignmentReader
        /// </summary>
        /// <param name="sr"></param>
        public PairwiseNucleotideSequenceAlignmentReader(StreamReader sr)
        {
            this.sr = sr;
            this.bufferSize = 5000;
        }

        /// <summary>
        /// Reads from a file containing Nucleotide sequence alignments
        /// </summary>
        /// <returns>A list full of SSRs</returns>
        public IPairwiseNucleotideSequenceAlignmentContainer GetNextAlignment()
        {
            //Instantiate the Working variables
            ReadMode readMode = ReadMode.Nonsense;
            SubAlignmentBuilder subBuild = new SubAlignmentBuilder();
            AlignmentContainerBuilder contBuilder = new AlignmentContainerBuilder();
            string databaseParent = "";
            string queryParent = "";
            int? lengthDatabaseParent = null;
            int? lengthQueryParent = null;
            int countLine = 0;
            bool? plusPlusStrand = null;
            bool invalid = true;

            IPairwiseNucleotideSequenceAlignmentContainer toRet = null;


            //The List of alignments which should be returned



            string line = "";


            //As long as the Buffer is not full and there are lines in the file proceed!
            while (invalid && (line = sr.ReadLine()) != null)
            {

                line = line.Trim();
                if (!(readMode == ReadMode.SubAlignment) && line == "" || line.StartsWith("#")) continue;


                if (line.StartsWith(">>"))
                {
                    contBuilder = new AlignmentContainerBuilder();
                    readMode = ReadMode.GenerallInfo;
                }
                else if (line.StartsWith(">"))
                {
                    subBuild = new SubAlignmentBuilder(databaseParent, queryParent, plusPlusStrand.Value, lengthDatabaseParent, lengthQueryParent, line);
                    readMode = ReadMode.SubAlignment;
                    countLine = 0;
                }
                else if (line.StartsWith("<<"))
                {
                    toRet=(IPairwiseNucleotideSequenceAlignmentContainer)contBuilder.GetPairwiseAlignmentContainer();
                    invalid = false;
                    readMode = ReadMode.Nonsense;
                    databaseParent = "";
                    queryParent = "";
                    lengthDatabaseParent = null;
                    lengthQueryParent = null;
                    plusPlusStrand = null;
                }
                else if (line.StartsWith("<"))
                {
                    contBuilder.AddSubalignment(subBuild.GetAlignment());
                    readMode = ReadMode.Nonsense;
                    countLine = 0;
                }

                else if (readMode == ReadMode.GenerallInfo)
                {
                    //dbName: FBgn0040060
                    //queName: ES6HDNH01BNV18 LEN=110 QL=1 QR=89
                    //subAlgn: 1
                    //strand: +/+
                    //dbLength: 
                    //queLength: 

                    if (line.StartsWith("dbName:")) databaseParent = line.Substring(8, line.Length - 8);
                    else if (line.StartsWith("queName:")) queryParent = line.Substring(9, line.Length - 9);
                    else if (line.StartsWith("dbLength:")) lengthDatabaseParent = Convert.ToInt32(line.Split(new char[] { ' ' })[1]);
                    else if (line.StartsWith("queLength:")) lengthQueryParent = Convert.ToInt32(line.Split(new char[] { ' ' })[1]);
                    else if (line.StartsWith("strand:"))
                    {
                        string temp = line.Split(new char[] { ' ' })[1];
                        temp = temp.Trim();

                        if (temp == "+/+") plusPlusStrand = true;
                        else if (temp == "+/-") plusPlusStrand = false;
                        else throw new InvalidDataException("this should not happen");
                    }

                }
                else if (readMode == ReadMode.SubAlignment)
                {
                    countLine++;
                    if (countLine % 4 == 1)
                    {
                        subBuild.AddDatabasesequence(line);
                    }
                    else if (countLine % 4 == 3)
                    {
                        subBuild.AddQuerysequence(line);
                    }
                }

            }



            return toRet;

        }

        public List<IPairwiseNucleotideSequenceAlignmentContainer> GetAllAlignments()
        {
            List<IPairwiseNucleotideSequenceAlignmentContainer> alignmentsList = new List<IPairwiseNucleotideSequenceAlignmentContainer>();

#if DEBUG
            System.Diagnostics.Stopwatch stopw = new System.Diagnostics.Stopwatch();
            stopw.Start();
#endif

            IPairwiseNucleotideSequenceAlignmentContainer temp;
            while ((temp = GetNextAlignment()) != null)
            {
                alignmentsList.Add(temp);
            }


#if DEBUG
            stopw.Stop();
            System.Diagnostics.Debug.WriteLine(String.Format("{0} PairwiseNucleotideSequenceAlignmentReader; All {1} Alignments read in {2} ms", DateTime.Now.TimeOfDay, alignmentsList.Count, stopw.ElapsedMilliseconds));
#endif
            return alignmentsList;

        }

        public List<IPairwiseNucleotideSequenceAlignmentContainer> GetAlignments(int number)
        {
            List<IPairwiseNucleotideSequenceAlignmentContainer> alignmentsList = new List<IPairwiseNucleotideSequenceAlignmentContainer>();

#if DEBUG
            System.Diagnostics.Stopwatch stopw = new System.Diagnostics.Stopwatch();
            stopw.Start();
#endif

            IPairwiseNucleotideSequenceAlignmentContainer temp;
            while (alignmentsList.Count<number && (temp = GetNextAlignment()) != null)
            {
                alignmentsList.Add(temp);
            }


#if DEBUG
            stopw.Stop();
            System.Diagnostics.Debug.WriteLine(String.Format("{0} PairwiseNucleotideSequenceAlignmentReader; A batch of {1} Alignments read in {2} ms", DateTime.Now.TimeOfDay, alignmentsList.Count, stopw.ElapsedMilliseconds));
#endif
            return alignmentsList;
        }



    }



    /// <summary>
    /// Writes the search results in a concise manner in the tab-delimited format
    /// TODO: definitely needs some additional work
    /// </summary>
    public class NucleotideSequenceAlignmentWriter_TabDelimited
    {
        private StreamWriter sw;
        private bool headerWriten = false;
        public NucleotideSequenceAlignmentWriter_TabDelimited(StreamWriter sw)
        {
            this.sw = sw;

        }
        public void WriteAlignments(List<PairwiseNucleotideSequenceAlignment> alignments)
        {

#if DEBUG
            System.Diagnostics.Stopwatch stopw = new System.Diagnostics.Stopwatch();
            stopw.Start();
#endif
            if (!headerWriten)
            {
                //              0               1           2       3                   4                5             6               7         8               9           10             11              12             13
                sw.WriteLine("database_name\tquery_name\tstrand\tlength\tstart_pos_database\tend_pos_database\tstart_pos_query\tend_pos_query\tscore\tsimilarity_woGaps\tsimilarity_wGaps\ttransitions\ttransversion\tsignificator"); //9
                headerWriten = true;
            }
            foreach (PairwiseNucleotideSequenceAlignment na in alignments)
            {
                sw.WriteLine(String.Format("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9:f}\t{10:f}\t{11}\t{12}\t{13}"
                    , na.DatabaseParent.Replace('\t', ' '),
                    na.QueryParent.Replace('\t', ' '),
                    na.PlusPlusStrand == true ? "+/+" : "+/-",
                    na.LengthAlignmentWithGaps,
                    na.StartDatabase,
                    na.EndDatabase,
                    na.StartQuery,
                    na.EndQuery,
                    na.Score,
                    na.SimilarityWithoutGaps,
                    na.SimilarityWithGaps,
                    na.Transitions,
                    na.Transversions,
                    na.Significator
                    ));
            }

#if DEBUG
            stopw.Stop();
            System.Diagnostics.Debug.WriteLine(String.Format("{0} NucleotideSequenceAlignmentWriter_TabDelimited; {1} Alignments written in TabDelimited style output format; Time used: {2} ms", DateTime.Now.TimeOfDay, alignments.Count, stopw.ElapsedMilliseconds));

#endif

        }



    }


   
  
}
