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


namespace Bio.Blast.Misc
{

    using Bio.Blast;
    using Bio.Seq;
    using Bio.Alignment;


    public interface IIntronBoundaryPolishing
    {
        System.Collections.Generic.List<Bio.Alignment.PairwiseNucleotideSequenceAlignment> GetPolishedAlignments(System.Collections.Generic.List<Bio.Alignment.PairwiseNucleotideSequenceAlignment> toPolish);
    }

    /// <summary>
    /// Refines the intron boundaries for a given collection of pairwise nucleotide sequence alignments;
    /// Note that the alignments have to have identical parents and plus plus strand
    /// </summary>
    public class IntronBoundaryPolishing_KeepBestOverlap : IIntronBoundaryPolishing
    {
        protected struct ExonPair
        {
            public PairwiseNucleotideSequenceAlignment pw3p;
            public PairwiseNucleotideSequenceAlignment pw5p;

        }
        private int minimumLengthOfPartial = 0;
        protected SubstitutionMatrix substitutionMatrix;
        System.IO.StreamWriter sw;

        public IntronBoundaryPolishing_KeepBestOverlap(SubstitutionMatrix substitutionMatrix)
        {
            this.substitutionMatrix = substitutionMatrix;


#if DEBUGh
                         sw= new System.IO.StreamWriter("G:\\Temp\\Report.txt", false, Encoding.ASCII);
#endif
        }

        public IntronBoundaryPolishing_KeepBestOverlap(SubstitutionMatrix substitutionMatrix, int minimumLengthOfPartialAlignments)
        {
            this.substitutionMatrix = substitutionMatrix;
            this.minimumLengthOfPartial = minimumLengthOfPartialAlignments;



        }

        //Destruktor!!
        ~IntronBoundaryPolishing_KeepBestOverlap()
        {
#if DEBUGh
                        if (sw != null)
                        {
                            try
                            {

                                sw.Close();
                            }
                            catch
                            {
                            }
                        }
#endif
        }

        public List<PairwiseNucleotideSequenceAlignment> GetPolishedAlignments(List<PairwiseNucleotideSequenceAlignment> toPolish)
        {
            //In case to polish contains 1 or less entries return it unaltered
            if (toPolish.Count < 2) return toPolish;

            foreach (PairwiseNucleotideSequenceAlignment pwa in toPolish)
            {
                if (pwa.DatabaseParent != toPolish[0].DatabaseParent || pwa.PlusPlusStrand != toPolish[0].PlusPlusStrand) throw new InvalidOperationException("All alignments have to have the same parents and must be from the same strand");
            }
            //Sort the list of possible exons new BlastUtility.HSPSort_DatabasePosition_QueryPosition<PairwiseNucleotideSequenceAlignment>()
            toPolish.Sort(new HSPSort_DatabasePosition_QueryPosition<PairwiseNucleotideSequenceAlignment>());

            List<PairwiseNucleotideSequenceAlignment> toRet = new List<PairwiseNucleotideSequenceAlignment>();

            //First fill a stack in a way that the 5' partial alignments will be first poped and the 3' will be last poped
            Stack<PairwiseNucleotideSequenceAlignment> stack = new Stack<PairwiseNucleotideSequenceAlignment>();
            for (int i = toPolish.Count - 1; i >= 0; i--)
            {
                stack.Push(toPolish[i]);
            }


            while (stack.Count > 1) //Pretty cool algorithm, 
            {
                ExonPair epToPolish = new ExonPair();
                epToPolish.pw5p = stack.Pop();
                epToPolish.pw3p = stack.Pop();

                ExonPair polished = Polish(epToPolish);

                if (polished.pw3p == null) //it is possible the the end polishing leaves short partial alginments empty, that is without any sequence, in that case use the remaining sequence and polish it with the next one in the stack
                {
                    stack.Push(polished.pw5p);
                }
                else if (polished.pw5p == null)
                {
                    stack.Push(polished.pw3p);
                }
                else
                {
                    stack.Push(polished.pw3p);
                    toRet.Add(polished.pw5p);
                }

            }
            if (stack.Count == 1) toRet.Add(stack.Pop());


            return toRet;
        }


        private ExonPair Polish(ExonPair toPolishPair)
        {
            PairwiseNucleotideSequenceAlignment leftAlignment = toPolishPair.pw5p;
            PairwiseNucleotideSequenceAlignment rightAlignment = toPolishPair.pw3p;

            if (leftAlignment.EndQuery < rightAlignment.StartQuery && leftAlignment.EndDatabase < rightAlignment.StartDatabase)
            {
                return toPolishPair;
            }
            // if (leftAlignment.EndDatabase > rightAlignment.StartDatabase) throw new Exception("should not happen");

            ///
            ///Create the working variables
            ///
            //0123456    45678901      
            //AAAATTT    TTTGGGGG
            int overlapQuery = leftAlignment.EndQuery - rightAlignment.StartQuery + 1;
            if (overlapQuery < 1) throw new Exception();

            ISequenceContainer rightquery = rightAlignment.Alignment.QuerySequence;
            ISequenceContainer rightDatabase = rightAlignment.Alignment.DatabaseSequence;
            ISequenceContainer leftQuery = leftAlignment.Alignment.QuerySequence;
            ISequenceContainer leftDatabase = leftAlignment.Alignment.DatabaseSequence;

            PairwiseAlignmentBuilder leftBuilder = new PairwiseAlignmentBuilder(leftAlignment.Alignment);
            PairwiseAlignmentBuilder rightBuilder = new PairwiseAlignmentBuilder(rightAlignment.Alignment);
            PairwiseAlignmentBuilder leftOverlapBuilder = new PairwiseAlignmentBuilder();
            PairwiseAlignmentBuilder rightOverlapBuilder = new PairwiseAlignmentBuilder();

            //Count the left and right overlap;
            int count_left_query = 0;
            int count_left_database = 0;
            int count_ToRemoveLeft = 0;
            int count_right_query = 0;
            int count_right_database = 0;
            int count_ToRemoveRight = 0;


            //Remove left
            for (int i = leftAlignment.Alignment.Length - 1; i >= 0; i--)
            {
                if (leftQuery[i] != '-') count_left_query++;
                if (leftDatabase[i] != '-') count_left_database++;

                leftOverlapBuilder.Append_5_prime(leftDatabase[i], leftQuery[i]);
                count_ToRemoveLeft++;
                if (count_left_query >= overlapQuery)
                {
                    if (i > 0)
                    {
                        if ((leftQuery[i - 1] == leftDatabase[i - 1]) //break only if the alingment does not end with a mismatch or a gap
                        && leftQuery[i - 1] != '-'
                        && leftDatabase[i - 1] != '-') break;
                    }
                    else break;

                }


            }
            leftBuilder.Remove_3_prime(count_ToRemoveLeft);



            //Remove right
            for (int i = 0; i < rightAlignment.Alignment.Length; i++)
            {
                if (rightquery[i] != '-') count_right_query++;
                if (rightDatabase[i] != '-') count_right_database++;

                rightOverlapBuilder.Append_3_prime(rightDatabase[i], rightquery[i]);
                count_ToRemoveRight++;
                if (count_right_query >= overlapQuery)
                {
                    if (i < rightAlignment.Alignment.Length - 1)
                    {
                        if ((rightquery[i + 1] == rightDatabase[i + 1])
                       && rightDatabase[i + 1] != '-'
                       && rightquery[i + 1] != '-') break;
                    }
                    else break;
                }
            }
            rightBuilder.Remove_5_prime(count_ToRemoveRight);


            PairwiseAlignment leftOverhangAlignment = leftOverlapBuilder.GetAlignment();
            PairwiseAlignment rightOverhangAlignment = rightOverlapBuilder.GetAlignment();
            float leftScore = substitutionMatrix.GetScoreForAlignment(leftOverhangAlignment);
            float rightScore = substitutionMatrix.GetScoreForAlignment(rightOverhangAlignment);

            ExonPair toRet = new ExonPair();
            if (leftScore >= rightScore)
            {
                toRet.pw5p = leftAlignment;

                PairwiseAlignment rightPwAl = rightBuilder.GetAlignment();
                PairwiseNucleotideSequenceAlignment rightAl = new PairwiseNucleotideSequenceAlignment(rightPwAl,
                 rightAlignment.DatabaseParent
                 , rightAlignment.QueryParent
                 , rightAlignment.StartDatabase + count_right_database
                 , rightAlignment.StartQuery + count_right_query
                 , rightAlignment.EndDatabase
                 , rightAlignment.EndQuery);
                rightAl.PlusPlusStrand = rightAlignment.PlusPlusStrand;
                rightAl.Score = rightAlignment.Score - rightScore;
                rightAl.SubstitutionMatrix = rightAlignment.SubstitutionMatrix;
                rightAl.AlgorithmUsedForAlignment = rightAlignment.AlgorithmUsedForAlignment;
                rightAl.LengthDatabaseParent = rightAlignment.LengthDatabaseParent;
                rightAl.LengthQueryParent = rightAlignment.LengthQueryParent;

                if (rightPwAl.DatabaseSequence.Length >= minimumLengthOfPartial) //Add the shortened sequence only if it is longer than zero 
                {
                    toRet.pw3p = rightAl;
                }
                else toRet.pw3p = null;

            }
            else
            {

                toRet.pw3p = rightAlignment;

                PairwiseAlignment leftPwAl = leftBuilder.GetAlignment();
                PairwiseNucleotideSequenceAlignment leftNsPwa = new PairwiseNucleotideSequenceAlignment(leftPwAl,
                 leftAlignment.DatabaseParent, leftAlignment.QueryParent,
                 leftAlignment.StartDatabase, leftAlignment.StartQuery,
                 leftAlignment.EndDatabase - count_left_database
                , leftAlignment.EndQuery - count_left_query);
                leftNsPwa.PlusPlusStrand = leftAlignment.PlusPlusStrand;
                leftNsPwa.Score = leftAlignment.Score - leftScore;
                leftNsPwa.SubstitutionMatrix = leftAlignment.SubstitutionMatrix;
                leftNsPwa.AlgorithmUsedForAlignment = leftAlignment.AlgorithmUsedForAlignment;
                leftNsPwa.LengthDatabaseParent = leftAlignment.LengthDatabaseParent;
                leftNsPwa.LengthQueryParent = leftAlignment.LengthQueryParent;

                if (leftPwAl.DatabaseSequence.Length >= minimumLengthOfPartial) //Add the shortened sequence only if it is longer than the minimum length
                {
                    toRet.pw5p = leftNsPwa;
                }
                else toRet.pw5p = null;
            }
#if DEBUGh
                        if (true)
                        {

                            sw.WriteLine("Left overlap");
                            sw.WriteLine(leftOverhangAlignment.DatabaseSequence.ToString());
                            sw.WriteLine(leftOverhangAlignment.SimilaritySequence.ToString());
                            sw.WriteLine(leftOverhangAlignment.QuerySequence.ToString());
                            sw.WriteLine();
                            sw.WriteLine("Right overlap");
                            sw.WriteLine(rightOverhangAlignment.DatabaseSequence.ToString());
                            sw.WriteLine(rightOverhangAlignment.SimilaritySequence.ToString());
                            sw.WriteLine(rightOverhangAlignment.QuerySequence.ToString());
                            sw.WriteLine();

                            if (leftScore >= rightScore)
                            {
                                sw.WriteLine("Kept left overlap");
                            }
                            else
                            {
                                sw.WriteLine("Kept right overlap");
                            }
                            sw.WriteLine();
                            sw.WriteLine();
                            sw.WriteLine();

   

                        }
#endif

            return toRet;



        }

    }



    public class IntronBoundaryPolishing_Cutoff : IIntronBoundaryPolishing
    {

        public IntronBoundaryPolishing_Cutoff()
        {
        }

        public List<PairwiseNucleotideSequenceAlignment> GetPolishedAlignments(List<PairwiseNucleotideSequenceAlignment> toPolish)
        {
            //pwal is a proven list of sequences which are on the same parent, have the same orientation
            //are in ascending order in the database sequence and in ascending order in the query sequence
            List<PairwiseNucleotideSequenceAlignment> temp = new List<PairwiseNucleotideSequenceAlignment>(toPolish.Count);
            toPolish.Sort(new Bio.Blast.Misc.HSPSort_DatabasePosition_QueryPosition<PairwiseNucleotideSequenceAlignment>());

            temp.Add(toPolish[0]);
            int lastQueryEnd = toPolish[0].EndQuery;
            for (int i = 1; i < toPolish.Count; i++)
            {

                if (lastQueryEnd >= toPolish[i].StartQuery)
                {
                    //If the sub alignment is greater than zero add it
                    PairwiseNucleotideSequenceAlignment pw = toPolish[i].SubAlignment_RelativeToQuery(lastQueryEnd + 1);
                    if (pw.Alignment.Length > 0) temp.Add(pw);

                }
                else
                {
                    temp.Add(toPolish[i]);
                }
                lastQueryEnd = toPolish[i].EndQuery;
            }
            return temp;

        }

    }

    public class BlastSignificator_808080 : ISignificator
    {
        private class SimilaritySlider
        {
            private int countSimilar = 0;

            public SimilaritySlider()
            {
            }

            public void Add(char c)
            {
                if (c == '|') countSimilar++;
            }

            public void Remove(char c)
            {
                if (c == '|') countSimilar--;
            }

            public float SimilarityPercent()
            {
                return (float)(countSimilar * 100) / 80.0F;
            }

        }
        public BlastSignificator_808080()
        {

        }


        public bool IsSignificant(IPairwiseAlignmentContainer simpleAlignment)
        {

            if (simpleAlignment.LengthAlignedDatabase < 80 || simpleAlignment.LengthAlignedQuery < 80) return false;


            ISequenceContainer similarity = simpleAlignment.Alignment.SimilaritySequence;
            SimilaritySlider slider = new SimilaritySlider();
            for (int i = 0; i < 80; i++)
            {
                slider.Add(similarity[i]);
            }
            if (slider.SimilarityPercent() >= 80.0F) return true;


            for (int i = 80; i < similarity.Length; i++)
            {
                slider.Add(similarity[i]);
                slider.Remove(similarity[i - 80]);
                if (slider.SimilarityPercent() >= 80.0F) return true;


            }
            return false;

        }

    }



        /// <summary>
        /// 1..according to the score
        /// </summary>
       internal class SortBlastSeed_Score : IComparer<BlastSeed>
        {
            public SortBlastSeed_Score()
            {
            }

            int IComparer<BlastSeed>.Compare(BlastSeed seed1, BlastSeed seed2)
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




            /// <summary>
            /// Sort PairwiseNucleotideSequenceAlignment; Descending according to the score
            /// </summary>
            public class HSPSort_ScoreDescending<TPairwise> : IComparer<TPairwise> where TPairwise : IPairwiseAlignmentContainer
            {
                public HSPSort_ScoreDescending()
                {
                }

                int IComparer<TPairwise>.Compare(TPairwise c1, TPairwise c2)
                {
                    //Sort in a threefold manner;
                    //First according to the databaseID;
                    //Second according to the databasePosition;
                    //and third according to the queryPosition;
                    if (c1.Score < c2.Score)
                    {
                        return -1;
                    }
                    else if (c1.Score > c2.Score)
                    {
                        return 1;
                    }
                    else
                    {
                        return 0;
                    }
                }
            }

            /// <summary>
            /// Sort PairwiseNucleotideSequenceAlignments; 
            /// 1..According to the databaseName
            /// 2..According to the score;
            /// </summary>
            public class HSPSort_DatabaseID_PlusPlus_DatabasePosition<TPairwise> : IComparer<TPairwise> where TPairwise : IPairwiseAlignmentContainer
            {
                public HSPSort_DatabaseID_PlusPlus_DatabasePosition()
                {
                }

                int IComparer<TPairwise>.Compare(TPairwise c1, TPairwise c2)
                {
                    //Sort in a threefold manner;
                    //First according to the databaseID;
                    //Second according to the databasePosition;
                    //and third according to the queryPosition;
                    int comp = String.Compare(c1.DatabaseParent, c2.DatabaseParent);

                    if (comp == 0)
                    {
                        if (c1.PlusPlusStrand == c2.PlusPlusStrand)
                        {
                            if (c1.StartDatabase < c2.StartDatabase)
                            {
                                return -1;
                            }
                            else if (c1.StartDatabase > c2.StartDatabase)
                            {
                                return 1;
                            }
                            else
                            {
                                return 0;
                            }
                        }
                        else if (c1.PlusPlusStrand == true) return -1;
                        else return 1;
                    }
                    else
                    {

                        return comp;

                    }
                }
            }




        }


        /// <summary>
        /// Sort Blast Seeds
        /// 1..according to the database id
        /// 1..according to the offset
        /// </summary>
        internal class SortBlastSeed_DatabaseID_Offset : IComparer<BlastSeed>
        {
            public SortBlastSeed_DatabaseID_Offset()
            {
            }

            int IComparer<BlastSeed>.Compare(BlastSeed seed1, BlastSeed seed2)
            {
                //Sort in a threefold manner;
                //First according to the databaseID;
                //Second according to the databasePosition;
                //and third according to the queryPosition;
                if (seed1.DatabaseID < seed2.DatabaseID)
                {
                    return -1;
                }
                else if (seed1.DatabaseID > seed2.DatabaseID)
                {
                    return 1;
                }
                else
                {
                    if (seed1.Offset < seed2.Offset) return -1;
                    else if (seed1.Offset > seed2.Offset) return 1;
                    else return 0;

                }
            }
        }

        /// <summary>
        /// Sort Blast Seeds
        /// 1..according to the database id
        /// 1..according to the database position
        /// </summary>
        internal class SortBlastSeed_DatabaseID_DatabasePosition : IComparer<BlastSeed>
        {
            public SortBlastSeed_DatabaseID_DatabasePosition()
            {
            }

            int IComparer<BlastSeed>.Compare(BlastSeed seed1, BlastSeed seed2)
            {
                //Sort in a threefold manner;
                //First according to the databaseID;
                //Second according to the databasePosition;
                if (seed1.DatabaseID < seed2.DatabaseID)
                {
                    return -1;
                }
                else if (seed1.DatabaseID > seed2.DatabaseID)
                {
                    return 1;
                }
                else
                {
                    if (seed1.DatabaseStartPosition < seed2.DatabaseStartPosition) return -1;
                    else if (seed1.DatabaseStartPosition > seed2.DatabaseStartPosition) return 1;
                    else return 0;

                }
            }
        }



        /// <summary>
        /// Sort Blast Seeds
        /// 1..according to the database id
        /// 2..according to the database position
        /// 3..according to the query position
        /// </summary>
        internal class SortBlastSeed_DatabaseID_DatabasePosition_QueryPosition : IComparer<BlastSeed>
        {
            public SortBlastSeed_DatabaseID_DatabasePosition_QueryPosition()
            {
            }

            int IComparer<BlastSeed>.Compare(BlastSeed seed1, BlastSeed seed2)
            {
                //Sort in a threefold manner;
                //First according to the databaseID;
                //Second according to the databasePosition;
                if (seed1.DatabaseID < seed2.DatabaseID)
                {
                    return -1;
                }
                else if (seed1.DatabaseID > seed2.DatabaseID)
                {
                    return 1;
                }
                else
                {
                    if (seed1.DatabaseStartPosition < seed2.DatabaseStartPosition) return -1;
                    else if (seed1.DatabaseStartPosition > seed2.DatabaseStartPosition) return 1;
                    else
                    {
                        if (seed1.QueryStartPosition < seed2.QueryStartPosition) return -1;
                        else if (seed1.QueryStartPosition > seed2.QueryStartPosition) return 1;
                        else return 0;
                    }

                }
            }
        }




        /// <summary>
        /// Sort PairwiseNucleotideSequenceAlignment; Descending according to the score
        /// </summary>
        internal class HSPSort_ScoreDescending<TPairwise> : IComparer<TPairwise> where TPairwise : IPairwiseAlignmentContainer
        {
            public HSPSort_ScoreDescending()
            {
            }

            int IComparer<TPairwise>.Compare(TPairwise c1, TPairwise c2)
            {
                //Sort in a threefold manner;
                //First according to the databaseID;
                //Second according to the databasePosition;
                //and third according to the queryPosition;
                if (c1.Score < c2.Score)
                {
                    return 1;
                }
                else if (c1.Score > c2.Score)
                {
                    return -1;
                }
                else
                {
                    return 0;
                }
            }
        }

        /// <summary>
        /// Sort PairwiseNucleotideSequenceAlignments; 
        /// 1..According to the databaseName
        /// 2..According to the score;
        /// </summary>
        internal class HSPSort_DatabaseID_PlusPlus_DatabasePosition_QueryPosition<TPairwise> : IComparer<TPairwise> where TPairwise : IPairwiseAlignmentContainer
        {
            public HSPSort_DatabaseID_PlusPlus_DatabasePosition_QueryPosition()
            {
            }

            int IComparer<TPairwise>.Compare(TPairwise c1, TPairwise c2)
            {
                //Sort in a threefold manner;
                //First according to the databaseID;
                //Second according to the databasePosition;
                //and third according to the queryPosition;
                int comp = String.Compare(c1.DatabaseParent, c2.DatabaseParent);

                if (comp == 0)
                {
                    if (c1.PlusPlusStrand == c2.PlusPlusStrand)
                    {
                        if (c1.StartDatabase < c2.StartDatabase)
                        {
                            return -1;
                        }
                        else if (c1.StartDatabase > c2.StartDatabase)
                        {
                            return 1;
                        }
                        else
                        {
                            if (c1.StartQuery < c2.StartQuery)
                            {
                                return -1;
                            }
                            else if (c1.StartQuery > c2.StartQuery)
                            {
                                return 1;
                            }
                            else return 0;


                        }
                    }
                    else if (c1.PlusPlusStrand == true) return -1;
                    else return 1;
                }
                else
                {

                    return comp;

                }
            }
        }


        /// <summary>
        /// 1..according to the database position
        /// 2..according to the query position
        /// </summary>
        internal class HSPSort_DatabasePosition_QueryPosition<TPairwise> : IComparer<TPairwise> where TPairwise : IPairwiseAlignmentContainer
        {
            public HSPSort_DatabasePosition_QueryPosition()
            {
            }

            int IComparer<TPairwise>.Compare(TPairwise hsp1, TPairwise hsp2)
            {

                if (hsp1.StartDatabase < hsp2.StartDatabase) return -1;
                else if (hsp1.StartDatabase > hsp2.StartDatabase) return 1;
                else
                {
                    if (hsp1.StartQuery < hsp2.StartQuery) return -1;
                    else if (hsp1.StartQuery > hsp2.StartQuery) return 1;
                    else return 0;
                }

            }

        }



        /// <summary>
        ///A seed of the BLAST implementation, actually is a positionable instance
        /// </summary>
        public class BlastSeed : IPositionable
        {
            /// <summary>
            /// Create a new Blast seed
            /// </summary>
            /// <param name="databasePos">the database start position</param>
            /// <param name="queryPos">the query start position</param>
            public BlastSeed(int databasePos, int queryPos)
            {
                this.databaseStartPos = databasePos;
                this.queryStartPos = queryPos;

            }
            public BlastSeed(int databaseID, int databasePos, int queryPos)
            {
                if (databaseID != 0) this.databaseID = (ushort)databaseID;
                this.databaseStartPos = databasePos;
                this.queryStartPos = queryPos;
            }


            public BlastSeed(int databaseID, int databasePos, int length, int queryPos)
            {
                if (databaseID != 0) this.databaseID = (ushort)databaseID;
                this.databaseStartPos = databasePos;
                this.databaseEndPos = databasePos + length - 1;
                this.queryStartPos = queryPos;
            }

            /// <summary>
            /// Create a new Blast Seed; Seeds are usually used to initiate BLAST alignments
            /// </summary>
            /// <param name="hashSeed">the BlastHashSeed which will be extended to contain the additionally required informations</param>
            /// <param name="queryPos">the startposition of the seed in the query sequence</param>
            /// <param name="length">the length of the BlastSeed, usually the wordsize</param>
            public BlastSeed(BlastHashSeed hashSeed, int queryPos)
            {
                this.queryStartPos = queryPos;
                if (hashSeed.DatabaseID != 0) this.databaseID = (ushort)hashSeed.DatabaseID;

                this.databaseStartPos = hashSeed.DatabasePosition;

            }
            //Fixed variables
            private ushort? databaseID = null;
            private int databaseStartPos;
            private int queryStartPos;

            //
            /// <summary>
            /// Flexible acquired variables, only used when absolutely needed
            /// </summary>
            private int? databaseEndPos = null;
            private int? queryEndPos = null;
            private int? offset = null;
            private float? score = null;


            /// <summary>
            /// Get or set the start position in the database sequence
            /// </summary>
            public int DatabaseStartPosition
            {
                get
                {
                    return this.databaseStartPos;
                }


            }



            /// <summary>
            /// Get the database id; this is an integer which can be used in an array of databases to retrieve the appropriate database sequence
            /// </summary>
            public int DatabaseID
            {
                get
                {
                    if (databaseID == null) return 0;
                    else return databaseID.Value;
                }
            }

            /// <summary>
            /// Get the start position of the query sequence
            /// </summary>
            public int QueryStartPosition
            {
                get
                {
                    return this.queryStartPos;
                }
            }



            public int DatabaseEndPosition
            {
                get
                {
                    return this.databaseEndPos.Value;
                }
                set
                {
                    this.databaseEndPos = value;
                }
            }

            public int QueryEndPosition
            {
                get
                {
                    return this.queryEndPos.Value;
                }
                set
                {
                    this.queryEndPos = value;
                }
            }
            /// <summary>
            /// Get or set the offset of a BLAST seed
            /// </summary>
            public int Offset
            {
                get
                {
                    return this.offset.Value;
                }
                set
                {
                    this.offset = value;
                }
            }

            /// <summary>
            /// Get or set the score of a Blast seed
            /// </summary>
            public float Score
            {
                get
                {
                    return this.score.Value;
                }
                set
                {
                    this.score = value;
                }
            }



            /// <summary>
            /// Get or set the start position of the BLAST seed
            /// </summary>
            public int? Start
            {
                get
                {
                    return this.DatabaseStartPosition;
                }
                set
                {
                    this.databaseStartPos = value.Value;
                }
            }

            /// <summary>
            /// Get or set the end position of the BLAST seed
            /// </summary>
            public int? End
            {
                get
                {
                    return this.databaseEndPos;
                }
                set
                {
                    this.databaseEndPos = value.Value;
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
                    throw new InvalidOperationException("A blast seed can never be the root");
                }
            }

            public string ParentName
            {
                get
                {
                    return this.DatabaseID.ToString();
                }
                set
                {
                    throw new Exception("The method or operation is not implemented.");
                }
            }


            public long Length
            {
                get
                {
                    return this.Start.Value - this.End.Value + 1;
                }
            }


            public string FeatureName
            {
                get { return "Blast Seed"; }
            }

        }


    public class BlastHashSeed
    {
        public BlastHashSeed(int databaseID, int databasePos)
        {
            this.databaseID = (ushort)databaseID;
            this.databasePos = databasePos;
        }
        public BlastHashSeed(int databasePos)
        {
            this.databasePos = databasePos;
        }

        private ushort? databaseID = null;
        private int databasePos;

        public int DatabaseID
        {
            get
            {
                if (databaseID == null) return 0;
                else return databaseID.Value;
            }
        }

        public int DatabasePosition
        {
            get
            {
                return databasePos;
            }
        }
    }



    /// <summary>
    /// Implementation of the SmithWatermanGotoh algorithm, one sided and overlappig version, should be used for large anchored alignments, like in a BLAST implementation;
    /// This means both Sequences start with an absolute high score at their respecitve 5'ends and are than extende in the 3' direction
    /// can be used when it is already known where the two sequences are similar, eg. BLAST implementations and subsequent 3'-extensions
    /// </summary>
    public class SmithWatermanGotoh_DynamicBanded_3p : SmithWatermanGotoh, IDynamicBandedDynamicProgramming
    {

        private bool highScoreAtTheBeginning = true;
        //private bool queryEnd = false;
        private int overlap;



        public SmithWatermanGotoh_DynamicBanded_3p(ISequenceContainer databaseSequence, ISequenceContainer querySequence, SubstitutionMatrix substitutionMatrix, int overhead)
        {
            this.penaltyGapExist = substitutionMatrix.GapExistPenalty;
            this.penaltyGapExtend = substitutionMatrix.GapExtendPenalty;
            this.databaseSequence = databaseSequence;
            this.querySequence = querySequence;
            this.substitutionMatrix = substitutionMatrix;
            this.overlap = overhead;


            InitializeMatrix();
        }

        /// <summary>
        /// Is the high score at the beginning of the alignment
        /// As long as the high score is not at the beginning of the alignment a better one might be possible
        /// </summary>
        public bool EndofDynamicExtension
        {
            get
            {
                return this.highScoreAtTheBeginning;
            }
        }






        /// <summary>
        /// Initializes the smith waterman alignment matrix, like in the original publication, testing for the highest Wk-insertions
        /// the highest Wk-deletions and than max(s(i,k),wk-del,wk-ins,0);
        /// </summary>
        protected override void InitializeMatrix()
        {
            //Initialize the matrix
            matrix = new MatrixPosition[databaseSequence.Length + 1, querySequence.Length + 1];

            //initialize the Gotoh arrays - containing the score of the best previous deletion
            //Max(Hik,Dk_1+u);
            this.Dk_1 = new float[databaseSequence.Length + 1];
            this.Qk_1 = new float[querySequence.Length + 1];

            //First the zeroPosition
            matrix[0, 0] = new MatrixPosition(0, LastMove.Diagonal);
            Dk_1[0] = 0;
            Qk_1[0] = 0;


            //Second set all values in the horizontal (database) to the initate value minus the penalties
            for (int i = 1; i < matrix.GetLength(0); i++)
            {
                float valueAttempt = 0 - penaltyGapExist - (penaltyGapExtend * (i - 1));

                matrix[i, 0] = new MatrixPosition(valueAttempt, LastMove.Deletion);
                //leafe as zero
                Dk_1[i] = valueAttempt - penaltyGapExist;
            }

            //Third set all values in the vertical (query) to zero
            for (int i = 1; i < matrix.GetLength(1); i++)
            {
                float valueAttempt = 0 - penaltyGapExist - (penaltyGapExtend * (i - 1));


                matrix[0, i] = new MatrixPosition(valueAttempt, LastMove.Insertion);
                Qk_1[i] = valueAttempt - penaltyGapExist;
            }

            //Go down use dimension of the query
            for (int k = 1; k < matrix.GetLength(1); k++)
            {
                for (int i = 1; i < matrix.GetLength(0); i++)
                {
                    //i=database sequence in the horizontal
                    //k=query sequence in the vertical


                    //the database sequence is in the horizontal, the query in the vertical axis of the matrix
                    //Diagonal score is the previous score and in addition the similarityValue;
                    float scoreDiagonal = matrix[i - 1, k - 1].score + substitutionMatrix.GetSimilarityValue(databaseSequence[i - 1], querySequence[k - 1]);

                    //Find the highest scoring insertion, testing all matrix to the upper side;
                    float downScoreInsertion = 0;
                    downScoreInsertion = Math.Max(matrix[i, k - 1].score - penaltyGapExist, Dk_1[i] - penaltyGapExtend);
                    Dk_1[i] = downScoreInsertion;


                    //Find the highest scroing deletion, testing all matrix entries to the left side
                    float rightScoreDeletion = 0;
                    rightScoreDeletion = Math.Max(matrix[i - 1, k].score - penaltyGapExist, Qk_1[k] - penaltyGapExtend);
                    Qk_1[k] = rightScoreDeletion;



                    matrix[i, k] = GetMaximumPosition(scoreDiagonal, downScoreInsertion, rightScoreDeletion);

#if DEBUG
                    MatrixPosition mp = GetMaximumPosition(scoreDiagonal, downScoreInsertion, rightScoreDeletion);
#endif

                    //Updating the highest scoring matrix entry
                    if (matrix[i, k].score > highScore)
                    {
                        //new highscore
                        highScore = matrix[i, k].score;
                        endAlignmentDatabase = i;
                        endAlignmentQuery = k;
                        this.highScoreAtTheBeginning = false;
                    }



                }
            }
            Traceback();
        }

        protected override void Traceback()
        {
            //If the matrix has not yet been constructed, initialize it
            if (matrix == null) InitializeMatrix();

            //Set the working variable, should the alignment be ignored until the overlap is reached or not
            bool ignoreAlignment = true;

            PairwiseAlignmentBuilder builder = new PairwiseAlignmentBuilder();
            //Actual position in the matrix; D..database; Q..query
            int posD = this.endAlignmentDatabase;
            int posQ = this.endAlignmentQuery;

            while (posD != 0 || posQ != 0)
            {

                //Move to the next character
                switch (matrix[posD, posQ].lastMove)
                {
                    case LastMove.Diagonal: //while
                        if (ignoreAlignment &&
                            ((databaseSequence.Length - posD) >= overlap) && ((querySequence.Length - posQ) >= overlap))
                        {
                            ignoreAlignment = false;
                            endAlignmentDatabase = posD;
                            endAlignmentQuery = posQ;
                            this.highScore = matrix[posD, posQ].score;
                        }
                        if (!ignoreAlignment) builder.Append_5_prime(databaseSequence[posD - 1], querySequence[posQ - 1]);
                        posD--; posQ--;

                        break;
                    case LastMove.Deletion:
                        int startPosD = posD;
                        do
                        {
                            if (!ignoreAlignment) builder.Append_5_prime(databaseSequence[posD - 1], '-');
                            posD--;

                        }
                        while (Math.Round(matrix[posD, posQ].score, 2) != Math.Round((matrix[startPosD, posQ].score + penaltyGapExist + penaltyGapExtend * (startPosD - posD - 1)), 2));

                        break;
                    case LastMove.Insertion:
                        int startPosQ = posQ;
                        do
                        {
                            if (!ignoreAlignment) builder.Append_5_prime('-', querySequence[posQ - 1]);
                            posQ--;

                        }
                        while (Math.Round(matrix[posD, posQ].score, 2) != Math.Round((matrix[posD, startPosQ].score + penaltyGapExist + penaltyGapExtend * (startPosQ - posQ - 1)), 2));


                        break;
                    default: throw new InvalidOperationException("this should not happen");
                }


            }


            if (endAlignmentDatabase > 0 && endAlignmentQuery > 0)
            {
                startAlignmentDatabase = posD + 1;
                startAlignmentQuery = posQ + 1;

                alignment = builder.GetAlignment();
            }

            //If the alignment length is zero, reset the values
            if (alignment != null && alignment.Length < 1)
            {
                alignment = null;
                this.highScoreAtTheBeginning = true;
                this.startAlignmentDatabase = 0;
                this.startAlignmentQuery = 0;
                this.endAlignmentDatabase = 0;
                this.endAlignmentQuery = 0;
            }
        }






    }


    public interface IDynamicBandedDynamicProgramming : IDynamicProgramming
    {
        bool EndofDynamicExtension { get;}

    }

    /// <summary>
    /// Implementation of the Smith Waterman algorithm, modified according to Gotoh
    /// This implementation (Gotoh) is much faster than the original Smith Waterman (M*M*N) and requires only a fraction of the memory (M*N)
    /// Returns only the highest scoring pair of sequences, the second and third best can not be obtained
    /// uses the equation: Wk= PExi - PExt * (x-1) for calculating the gap penalties
    /// </summary>
    public class SmithWatermanGotoh_DynamicBanded_5p : SmithWatermanGotoh, IDynamicBandedDynamicProgramming
    {
        private bool highScoreAtTheBeginning = true;
        //private bool queryEnd = false;
        private int overlap;



        public SmithWatermanGotoh_DynamicBanded_5p(ISequenceContainer databaseSequence, ISequenceContainer querySequence, SubstitutionMatrix substitutionMatrix, int overhead)
        {
            if (databaseSequence.Length < 1 || querySequence.Length < 1) throw new InvalidSequenceException("Error in Smith-Waterman-Gotoh dynamic programming. Length of database sequence or query sequence has to be at least 1 bp");


            this.penaltyGapExist = substitutionMatrix.GapExistPenalty;
            this.penaltyGapExtend = substitutionMatrix.GapExtendPenalty;
            this.databaseSequence = databaseSequence;
            this.querySequence = querySequence;
            this.substitutionMatrix = substitutionMatrix;
            this.overlap = overhead;


            InitializeMatrix();
        }



        /// <summary>
        /// Initializes the smith waterman alignment matrix, like in the original publication, testing for the highest Wk-insertions
        /// the highest Wk-deletions and than max(s(i,k),wk-del,wk-ins,0);
        /// </summary>
        protected override void InitializeMatrix()
        {
            //Initialize the matrix
            matrix = new MatrixPosition[databaseSequence.Length + 1, querySequence.Length + 1];

            //initialize the Gotoh arrays - containing the score of the best previous deletion
            //Max(Hik,Dk_1+u);
            this.Dk_1 = new float[databaseSequence.Length + 1];
            this.Qk_1 = new float[querySequence.Length + 1];



            //First set all values in the horizontal (database) to zero
            for (int i = 0; i < matrix.GetLength(0); i++)
            {
                matrix[i, 0] = new MatrixPosition(0, LastMove.None);
                Dk_1[i] = 0 - penaltyGapExist;
            }

            //Second set all values in the vertical (query) to zero
            for (int i = 0; i < matrix.GetLength(1); i++)
            {
                matrix[0, i] = new MatrixPosition(0, LastMove.None);
                Qk_1[i] = 0 - penaltyGapExist;
            }

            //Go down use dimension of the query
            for (int k = 1; k < matrix.GetLength(1); k++)
            {
                for (int i = 1; i < matrix.GetLength(0); i++)
                {
                    //i=database sequence in the horizontal
                    //k=query sequence in the vertical


                    //the database sequence is in the horizontal, the query in the vertical axis of the matrix
                    //Diagonal score is the previous score and in addition the similarityValue;
                    float scoreDiagonal = matrix[i - 1, k - 1].score + substitutionMatrix.GetSimilarityValue(databaseSequence[i - 1], querySequence[k - 1]);

                    //Find the highest scoring insertion, testing all matrix to the upper side;
                    float downScoreInsertion = 0;
                    downScoreInsertion = this.Max(matrix[i, k - 1].score - penaltyGapExist, Dk_1[i] - penaltyGapExtend);
                    Dk_1[i] = downScoreInsertion;


                    //Find the highest scroing deletion, testing all matrix entries to the left side
                    float rightScoreDeletion = 0;
                    rightScoreDeletion = Math.Max(matrix[i - 1, k].score - penaltyGapExist, Qk_1[k] - penaltyGapExtend);
                    Qk_1[k] = rightScoreDeletion;



                    matrix[i, k] = GetMaximumPosition(scoreDiagonal, downScoreInsertion, rightScoreDeletion, 0);

#if DEBUG
                    MatrixPosition mp = GetMaximumPosition(scoreDiagonal, downScoreInsertion, rightScoreDeletion, 0);
#endif

                    /*/Updating the highest scoring matrix entry
                                if (matrix[i, k].score > highScore)
                                {
                                    //new highscore
                                    highScore = matrix[i, k].score;
                                    endAlignmentDatabase = i;
                                    endAlignmentQuery = k;
                                }
                                 */



                }
            }
            //Call the traceback method
            Traceback();
        }


        protected override void Traceback()
        {

            PairwiseAlignmentBuilder builder = new PairwiseAlignmentBuilder();
            //Actual position in the matrix; D..database; Q..query
            int posD = this.databaseSequence.Length;
            int posQ = this.querySequence.Length;
            int lastPosD = posD;
            int lastPosQ = posQ;

            while (matrix[posD, posQ].lastMove != LastMove.None)
            {

                //Move to the next character
                switch (matrix[posD, posQ].lastMove)
                {
                    case LastMove.Diagonal: //while

                        posD--; posQ--;

                        if (posD >= overlap && posQ >= overlap)
                        {
                            lastPosD = posD;
                            lastPosQ = posQ;
                            builder.Append_5_prime(databaseSequence[posD], querySequence[posQ]);
                        }
                        else
                        {
                            goto endOfLoop;
                        }

                        break;
                    case LastMove.Deletion:
                        int startPosD = posD;
                        do
                        {
                            builder.Append_5_prime(databaseSequence[posD - 1], '-');
                            posD--;
                            if (!(posD > overlap))
                            {
                                builder.Remove_5_prime(startPosD - posD);
                                goto endOfLoop;
                            }
#if DEBUG
                            double a = Math.Round(matrix[posD, posQ].score, 4);
                            double b = Math.Round((matrix[startPosD, posQ].score + penaltyGapExist + penaltyGapExtend * (startPosD - posD - 1)), 4);
#endif

                        }
                        while (Math.Round(matrix[posD, posQ].score, 4) != Math.Round((matrix[startPosD, posQ].score + penaltyGapExist + penaltyGapExtend * (startPosD - posD - 1)), 4));

                        break;
                    case LastMove.Insertion:
                        int startPosQ = posQ;
                        do
                        {
                            builder.Append_5_prime('-', querySequence[posQ - 1]);
                            posQ--;
                            if (!(posQ > overlap))
                            {
                                builder.Remove_5_prime(startPosQ - posQ);
                                goto endOfLoop;
                            }
                        }
                        while (Math.Round(matrix[posD, posQ].score, 2) != Math.Round((matrix[posD, startPosQ].score + penaltyGapExist + penaltyGapExtend * (startPosQ - posQ - 1)), 2));

                        break;
                    default: throw new InvalidOperationException("this should not happen");
                }


            }

        endOfLoop: ;


            if (matrix[databaseSequence.Length, querySequence.Length].score != 0)
            {
                startAlignmentDatabase = -(databaseSequence.Length - lastPosD);
                startAlignmentQuery = -(querySequence.Length - lastPosQ);
                endAlignmentDatabase = -1;
                endAlignmentQuery = -1;
                this.highScore = matrix[databaseSequence.Length, querySequence.Length].score - matrix[lastPosD, lastPosQ].score;
                this.highScoreAtTheBeginning = false;
                alignment = builder.GetAlignment();
            }

            if (alignment != null && alignment.Length < 1)
            {
                alignment = null;
                this.highScoreAtTheBeginning = true;
                this.startAlignmentDatabase = 0;
                this.startAlignmentQuery = 0;
                this.endAlignmentDatabase = 0;
                this.endAlignmentQuery = 0;
            }
 


        }





        public override int Start_Database
        {
            get
            {
                //012345678901
                //AAAATTATATAT leng=12
                //123456789012
                //  | pos=3
                return this.startAlignmentDatabase;
            }
        }

        /// <summary>
        /// Start postion of the optimal alignment with respect to the query sequence
        /// </summary>
        public override int Start_Query
        {
            get
            {
                return this.startAlignmentQuery;
            }
        }

        /// <summary>
        /// End position of the optimal alignment with respect to the database sequence
        /// </summary>
        public override int End_Database
        {
            get
            {
                return this.endAlignmentDatabase;
            }
        }

        /// <summary>
        /// End position of the optimal alignment with respect to the query sequence
        /// </summary>
        public override int End_Query
        {
            get
            {
                return endAlignmentQuery;
            }
        }




        public bool EndofDynamicExtension
        {
            get
            {
                return this.highScoreAtTheBeginning;
            }
        }


    }


    /// <summary>
    /// Implementation of the Smith Waterman algorithm, modified according to Gotoh.
    /// Specialy adapted for the 454 sequencing, uses boundarie cross penaltys for gaps extending beyond poly-N tracts
    /// The direction of the 454 read is important, the query sequence has to have the 5'->3' direction 
    /// P..Plus/Plus
    /// </summary>
    public class SmithWatermanGotoh_DynamicBanded_454P_3p : SmithWatermanGotoh, IDynamicBandedDynamicProgramming
    {
        private float[] geD; //gap exist database
        private float[] geQ; //gap exist query
        private bool[] boundaryD;
        private bool[] boundaryQ;
        private float[] Dk_1_var;
        private bool[] dk_1_var_boundary;
        private bool[] qk_1_var_boundary;
        private float[] Qk_1_var;
        private float boundaryCrossPenalty;

        private bool endOfAlignment = true;

        private int overlap;



        protected SmithWatermanGotoh_DynamicBanded_454P_3p()
        {
        }

        public SmithWatermanGotoh_DynamicBanded_454P_3p(ISequenceContainer databaseSequence, ISequenceContainer querySequence, SubstitutionMatrix substitutionMatrix, int overhead)
        {
            if (databaseSequence.Length < 1 || querySequence.Length < 1) throw new InvalidSequenceException("Error in Smith-Waterman-Gotoh dynamic programming. Length of database sequence or query sequence has to be at least 1 bp");


            this.penaltyGapExist = substitutionMatrix.GapExistPenalty;
            this.penaltyGapExtend = substitutionMatrix.GapExtendPenalty;
            this.databaseSequence = databaseSequence;
            this.querySequence = querySequence;
            this.substitutionMatrix = substitutionMatrix;
            this.boundaryCrossPenalty = 2 * penaltyGapExtend;
            this.overlap = overhead;

            InitializeMatrix();
        }

        public SmithWatermanGotoh_DynamicBanded_454P_3p(ISequenceContainer databaseSequence, ISequenceContainer querySequence, SubstitutionMatrix substitutionMatrix, float homopolymereTransgressionPenalty, int overhead)
        {
            if (databaseSequence.Length < 1 || querySequence.Length < 1) throw new InvalidSequenceException("Error in Smith-Waterman-Gotoh dynamic programming. Length of database sequence or query sequence has to be at least 1 bp");


            this.penaltyGapExist = substitutionMatrix.GapExistPenalty;
            this.penaltyGapExtend = substitutionMatrix.GapExtendPenalty;
            this.databaseSequence = databaseSequence;
            this.querySequence = querySequence;
            this.substitutionMatrix = substitutionMatrix;
            this.boundaryCrossPenalty = homopolymereTransgressionPenalty;
            this.overlap = overhead;

            InitializeMatrix();
        }


        /// <summary>
        /// Initializes the smith waterman alignment matrix, like in the original publication, testing for the highest Wk-insertions
        /// the highest Wk-deletions and than max(s(i,k),wk-del,wk-ins,0);
        /// </summary>
        protected override void InitializeMatrix()
        {
            //Initialize the matrix
            matrix = new MatrixPosition[databaseSequence.Length + 1, querySequence.Length + 1];

            //initialize the Gotoh arrays - containing the score of the best previous deletion
            //Max(Hik,Dk_1+u);
            this.Dk_1 = new float[databaseSequence.Length + 1];
            this.Dk_1_var = new float[databaseSequence.Length + 1];
            this.dk_1_var_boundary = new bool[databaseSequence.Length + 1];
            this.qk_1_var_boundary = new bool[querySequence.Length + 1];
            this.Qk_1 = new float[querySequence.Length + 1];
            this.Qk_1_var = new float[querySequence.Length + 1];



            matrix[0, 0] = new MatrixPosition(0, LastMove.Diagonal);

            //First set all values in the horizontal (database) to zero
            for (int i = 1; i < matrix.GetLength(0); i++)
            {
                float valueAttempt = 0.0F - penaltyGapExist - (penaltyGapExtend * (i - 1));

                matrix[i, 0] = new MatrixPosition(valueAttempt, LastMove.Deletion);
                Dk_1[i] = valueAttempt - penaltyGapExist;
                Dk_1_var[i] = valueAttempt - penaltyGapExist;
            }


            //Second set all values in the vertical (query) to zero
            for (int i = 1; i < matrix.GetLength(1); i++)
            {
                float valueAttempt = 0.0F - penaltyGapExist - (penaltyGapExtend * (i - 1));
                matrix[0, i] = new MatrixPosition(valueAttempt, LastMove.Insertion);
                Qk_1[i] = valueAttempt - penaltyGapExist;
                Qk_1_var[i] = valueAttempt - penaltyGapExist;
            }

            /*
         for (int i = 1; i < matrix.GetLength(0); i++)
        {
            float valueAttempt = 0.0F - penaltyGapExist - (penaltyGapExtend * (i - 1));

            matrix[i, 0] = new MatrixPosition(valueAttempt, LastMove.Deletion);
            //leafe as zero
            Dk_1[i] = valueAttempt-penaltyGapExist;
        }

        //Third set all values in the vertical (query) to zero
        for (int i = 1; i < matrix.GetLength(1); i++)
        {
            float valueAttempt = 0.0F - penaltyGapExist - (penaltyGapExtend * (i - 1));


            matrix[0, i] = new MatrixPosition(valueAttempt, LastMove.Insertion);
            Qk_1[i] = valueAttempt-penaltyGapExist;
        }
             */

            ProcessGapExistPenalty(); //Create position specific gap existence penalty matrix

            //Go down use dimension of the query
            for (int k = 1; k < matrix.GetLength(1); k++)
            {
                for (int i = 1; i < matrix.GetLength(0); i++)
                {
                    float gp = Math.Min(geD[i], geQ[k]);

                    float scoreDiagonal = matrix[i - 1, k - 1].score + substitutionMatrix.GetSimilarityValue(databaseSequence[i - 1], querySequence[k - 1]);


                    //Find the highest scoring insertion, testing all matrix to the upper side;
                    Dk_1[i] = this.Max(matrix[i, k - 1].score - penaltyGapExist, Dk_1[i] - penaltyGapExtend);
                    Dk_1_var[i] = this.Max(matrix[i, k - 1].score - penaltyGapExtend, Dk_1_var[i] - penaltyGapExtend);
                    if (boundaryQ[k - 1] == true)
                    {
                        Dk_1_var[i] = Dk_1_var[i] - boundaryCrossPenalty;

                    }

                    float downScoreInsertion = Math.Max(Dk_1[i], (Dk_1_var[i] + penaltyGapExtend - gp));






                    //Find the highest scroing deletion, testing all matrix entries to the left side
                    Qk_1[k] = this.Max(matrix[i - 1, k].score - penaltyGapExist, Qk_1[k] - penaltyGapExtend);
                    Qk_1_var[k] = this.Max(matrix[i - 1, k].score - penaltyGapExtend, Qk_1_var[k] - penaltyGapExtend);
                    if (boundaryD[i - 1] == true)
                    {
                        Qk_1_var[k] = Qk_1_var[k] - boundaryCrossPenalty;

                    }

                    float rightScoreDeletion = Math.Max(Qk_1[k], (Qk_1_var[k] + penaltyGapExtend - gp));








                    matrix[i, k] = GetMaximumPosition(scoreDiagonal, downScoreInsertion, rightScoreDeletion);

#if DEBUG
                    MatrixPosition mp = GetMaximumPosition(scoreDiagonal, downScoreInsertion, rightScoreDeletion);
#endif

                    //Updating the highest scoring matrix entry
                    if (matrix[i, k].score > highScore)
                    {
                        //new highscore
                        highScore = matrix[i, k].score;
                        endAlignmentDatabase = i;
                        endAlignmentQuery = k;
                        this.endOfAlignment = false;
                    }



                }
            }
            //Call the traceback method
            Traceback();
        }

        private void ProcessGapExistPenalty()
        {

            float maxStep = Math.Abs(substitutionMatrix.HighestScore) + Math.Abs(substitutionMatrix.LowestScore) - 0.1F;
            float lowestGapExist = penaltyGapExist - maxStep;
            //if the lowestGapExist penalty is smaller than the gapExtend penalty increase it to the gapExtend penalty
            lowestGapExist = (lowestGapExist < penaltyGapExtend ? penaltyGapExtend + 0.1F : lowestGapExist);

            //Initialize the new arrays
            geD = new float[databaseSequence.Length + 1];
            geQ = new float[querySequence.Length + 1];
            boundaryD = new bool[databaseSequence.Length + 1];
            boundaryQ = new bool[querySequence.Length + 1];
            for (int i = 0; i < boundaryD.Length; i++)
            {
                boundaryD[i] = false;
            }
            for (int i = 0; i < boundaryQ.Length; i++)
            {
                boundaryQ[i] = false;
            }




            //The first two characters for every position specific gap penalty matrix are identical to the default gap existence penalty.
            geD[0] = this.penaltyGapExist;
            geQ[0] = this.penaltyGapExist;
            if (databaseSequence.Length >= 1) geD[1] = this.penaltyGapExist;
            if (querySequence.Length >= 1) geQ[1] = this.penaltyGapExist;

            //|0123456789
            //|TTTTTTTTTTTTT
            //99876543211111 geD
            //0123456789

            for (int i = 0; i < databaseSequence.Length - 1; i++)
            {

                if (databaseSequence[i] != databaseSequence[i + 1])
                {
                    geD[i + 2] = penaltyGapExist;
                    boundaryD[i + 1] = true;

                }
                else
                {
                    int k = 1; //k==1 -> unique k==2 AA; k==4 AAAA
                    while (i+k < databaseSequence.Length && databaseSequence[i] == databaseSequence[i + k])
                    {
                        k++;
                    }
                    float step = (penaltyGapExist - lowestGapExist) / (float)(k - 1);

                    //Fill until the new match starts
                    int l = 0;
                    for (; l < k - 1; l++)
                    {
                        geD[i + 2 + l] = geD[i + 1 + l] - step;
#if DEBUG
                        float assign = (geD[i + 1 + l] - step);
                        if (assign <= penaltyGapExtend) throw new Exception("Fucking impossible");
#endif
                    }
                    i = i + l - 1;


                }

            }


            //|0123456789
            //|TTTTTTTTTTTTT
            //99876543211111 geD
            //0123456789
            for (int i = 0; i < querySequence.Length - 1; i++)
            {

                if (querySequence[i] != querySequence[i + 1])
                {
                    geQ[i + 2] = penaltyGapExist;
                    boundaryQ[i + 1] = true;

                }
                else
                {
                    int k = 1; //k==1 -> unique k==2 AA; k==4 AAAA
                    while (i + k < querySequence.Length && querySequence[i] == querySequence[i + k])
                    {
                        k++;
                    }
                    float step = (penaltyGapExist - lowestGapExist) / (float)(k - 1);

                    //Fill until the new match starts
                    int l = 0;
                    for (; l < k - 1; l++)
                    {
                        geQ[i + 2 + l] = geQ[i + 1 + l] - step;

#if DEBUG
                        if ((geQ[i + 1 + l] - step) <= penaltyGapExtend) throw new Exception("Fucking impossible");
#endif
                    }

                    i = i + l - 1;


                }

            }





        }


        protected override void Traceback()
        {

            PairwiseAlignmentBuilder builder = new PairwiseAlignmentBuilder();
            //Actual position in the matrix; D..database; Q..query
            int posD = this.endAlignmentDatabase;
            int posQ = this.endAlignmentQuery;
            float gp = penaltyGapExist;
            bool predecesingGap = false;

            //As long as the required overlap has not been reached the alignment has to be ignored
            bool ignoreAlignment = true;

            while (posD != 0 && posQ != 0)
            {
                switch (matrix[posD, posQ].lastMove)
                {
                    case LastMove.Diagonal: //while
                        if (ignoreAlignment &&
                         ((databaseSequence.Length - posD) >= overlap) && ((querySequence.Length - posQ) >= overlap) && !predecesingGap)// && ((posD==0||posQ==0)||( boundaryD[posD - 1] == true || boundaryQ[posQ - 1] == true))
                        {
                            ignoreAlignment = false;
                            endAlignmentDatabase = posD;
                            endAlignmentQuery = posQ;
                            this.highScore = matrix[posD, posQ].score;
                        }
                        if (!ignoreAlignment) builder.Append_5_prime(databaseSequence[posD - 1], querySequence[posQ - 1]);
                        posD--; posQ--;
                        predecesingGap = false;
                        break;
                    case LastMove.Deletion:
                        int startPosD = posD;
                        gp = Math.Min(geD[posD], geQ[posQ]);
                        do
                        {
                            if (!ignoreAlignment) builder.Append_5_prime(databaseSequence[posD - 1], '-');
                            if (gp != penaltyGapExist && boundaryD[posD - 1] == true)
                            {
                                gp = gp + boundaryCrossPenalty;
                                if (gp > penaltyGapExist) gp = penaltyGapExist;
                            }
                            posD--;
#if DEBUG
                            double a = Math.Round(matrix[posD, posQ].score, 4);
                            double b = Math.Round((matrix[startPosD, posQ].score + gp + penaltyGapExtend * (startPosD - posD - 1)), 4);
#endif
                            predecesingGap = true;
                        }
                        while (Math.Abs(Math.Round(matrix[posD, posQ].score, 6) - Math.Round((matrix[startPosD, posQ].score + gp + penaltyGapExtend * (startPosD - posD - 1)), 6)) > 0.0001);
                        break;
                    case LastMove.Insertion:
                        int startPosQ = posQ;
                        gp = Math.Min(geD[posD], geQ[posQ]);
                        do
                        {
                            if (!ignoreAlignment) builder.Append_5_prime('-', querySequence[posQ - 1]);
                            if (gp != penaltyGapExist && boundaryQ[posQ - 1] == true)
                            {
                                gp = gp + boundaryCrossPenalty;
                                if (gp > penaltyGapExist) gp = penaltyGapExist;
                            }
                            posQ--;
#if DEBUG
                            double a = Math.Round(matrix[posD, posQ].score, 4);
                            double b = Math.Round((matrix[posD, startPosQ].score + gp + penaltyGapExtend * (startPosQ - posQ - 1)), 4);
#endif
                            predecesingGap = true;
                        }
                        while (Math.Abs(Math.Round(matrix[posD, posQ].score, 6) - Math.Round((matrix[posD, startPosQ].score + gp + penaltyGapExtend * (startPosQ - posQ - 1)), 6)) > 0.0001);
                        break;
                    default: throw new InvalidOperationException("this should not happen");
                }
            }

            if (endAlignmentDatabase > 0 && endAlignmentQuery > 0)
            {
                startAlignmentDatabase = posD + 1;
                startAlignmentQuery = posQ + 1;
                alignment = builder.GetAlignment();
   
            }
            //If the alignment length is zero, reset the values
            if (alignment!=null && alignment.Length < 1)
            {
                alignment = null;
                this.endOfAlignment = true;
                this.startAlignmentDatabase = 0;
                this.startAlignmentQuery = 0;
                this.endAlignmentDatabase = 0;
                this.endAlignmentQuery = 0;
            }
        }





        public string ToStringVariableGapPenaltyDatabase()
        {
            StringBuilder sb = new StringBuilder("");
            for (int i = 0; i < geD.Length; i++)
            {
                sb.Append((int)geD[i]);
            }
            return sb.ToString();
        }

        public string ToStringVariableGapPenaltyQuery()
        {

            StringBuilder sb = new StringBuilder("");
            for (int i = 0; i < geQ.Length; i++)
            {
                sb.Append((int)geQ[i]);
            }
            return sb.ToString();

        }

        public string ToStringBoundariesDatabase()
        {

            StringBuilder sb = new StringBuilder("");
            for (int i = 0; i < boundaryD.Length; i++)
            {
                sb.Append(boundaryD[i] == true ? "T" : "F");
            }
            return sb.ToString();

        }

        public string ToStringBoundariesQuery()
        {

            StringBuilder sb = new StringBuilder("");
            for (int i = 0; i < boundaryQ.Length; i++)
            {
                sb.Append(boundaryQ[i] == true ? "T" : "F");
            }
            return sb.ToString();

        }




        #region IDynamicBandedDynamicProgramming Members

        public bool EndofDynamicExtension
        {
            get { return this.endOfAlignment; }
        }

        #endregion
    }


    /// <summary>
    /// Implementation of the Smith Waterman algorithm, modified according to Gotoh.
    /// Specialy adapted for the 454 sequencing, uses boundarie cross penaltys for gaps extending beyond poly-N tracts
    /// The direction of the 454 read is important, the query sequence has to have the 5'->3' direction 
    /// P..Plus/Plus
    /// </summary>
    public class SmithWatermanGotoh_DynamicBanded_454P_5p : SmithWatermanGotoh, IDynamicBandedDynamicProgramming
    {
        private float[] geD; //gap exist database
        private float[] geQ; //gap exist query
        private bool[] boundaryD;
        private bool[] boundaryQ;
        private int overlap;

        private float[] Dk_1_var;
        private bool[] dk_1_var_boundary;
        private bool[] qk_1_var_boundary;
        private float[] Qk_1_var;
        private float gapLength = 3.0F;
        private float boundaryCrossPenalty;

        private bool highScoreAtTheBeginning = true;



        protected SmithWatermanGotoh_DynamicBanded_454P_5p()
        {
        }

        public SmithWatermanGotoh_DynamicBanded_454P_5p(ISequenceContainer databaseSequence, ISequenceContainer querySequence, SubstitutionMatrix substitutionMatrix, int overhead)
        {
            if (databaseSequence.Length < 1 || querySequence.Length < 1) throw new InvalidSequenceException("Error in Smith-Waterman-Gotoh dynamic programming. Length of database sequence or query sequence has to be at least 1 bp");


            this.penaltyGapExist = substitutionMatrix.GapExistPenalty;
            this.penaltyGapExtend = substitutionMatrix.GapExtendPenalty;
            this.databaseSequence = databaseSequence;
            this.querySequence = querySequence;
            this.substitutionMatrix = substitutionMatrix;
            this.boundaryCrossPenalty = 2 * penaltyGapExtend;
            this.overlap = overhead;

            InitializeMatrix();
        }

        public SmithWatermanGotoh_DynamicBanded_454P_5p(ISequenceContainer databaseSequence, ISequenceContainer querySequence, SubstitutionMatrix substitutionMatrix, float homopolymereTransgressionPenalty, int overhead)
        {
            if (databaseSequence.Length < 1 || querySequence.Length < 1) throw new InvalidSequenceException("Error in Smith-Waterman-Gotoh dynamic programming. Length of database sequence or query sequence has to be at least 1 bp");


            this.penaltyGapExist = substitutionMatrix.GapExistPenalty;
            this.penaltyGapExtend = substitutionMatrix.GapExtendPenalty;
            this.databaseSequence = databaseSequence;
            this.querySequence = querySequence;
            this.substitutionMatrix = substitutionMatrix;
            this.boundaryCrossPenalty = homopolymereTransgressionPenalty;
            this.overlap = overhead;

            InitializeMatrix();
        }


        /// <summary>
        /// Initializes the smith waterman alignment matrix, like in the original publication, testing for the highest Wk-insertions
        /// the highest Wk-deletions and than max(s(i,k),wk-del,wk-ins,0);
        /// </summary>
        protected override void InitializeMatrix()
        {
            //Initialize the matrix
            matrix = new MatrixPosition[databaseSequence.Length + 1, querySequence.Length + 1];

            //initialize the Gotoh arrays - containing the score of the best previous deletion
            //Max(Hik,Dk_1+u);
            this.Dk_1 = new float[databaseSequence.Length + 1];
            this.Dk_1_var = new float[databaseSequence.Length + 1];
            this.dk_1_var_boundary = new bool[databaseSequence.Length + 1];
            this.qk_1_var_boundary = new bool[querySequence.Length + 1];
            this.Qk_1 = new float[querySequence.Length + 1];
            this.Qk_1_var = new float[querySequence.Length + 1];



            //First set all values in the horizontal (database) to zero
            for (int i = 0; i < matrix.GetLength(0); i++)
            {
                matrix[i, 0] = new MatrixPosition(0, LastMove.None);
                Dk_1[i] = 0 - penaltyGapExist;
                Dk_1_var[i] = 0;
            }


            //Second set all values in the vertical (query) to zero
            for (int i = 0; i < matrix.GetLength(1); i++)
            {
                matrix[0, i] = new MatrixPosition(0, LastMove.None);
                Qk_1[i] = 0 - penaltyGapExist;
                Qk_1_var[i] = 0;
            }

            ProcessGapExistPenalty(); //Create position specific gap existence penalty matrix

            //Go down use dimension of the query
            for (int k = 1; k < matrix.GetLength(1); k++)
            {
                for (int i = 1; i < matrix.GetLength(0); i++)
                {
                    float gp = Math.Min(geD[i], geQ[k]);

                    float scoreDiagonal = matrix[i - 1, k - 1].score + substitutionMatrix.GetSimilarityValue(databaseSequence[i - 1], querySequence[k - 1]);


                    //Find the highest scoring insertion, testing all matrix to the upper side;
                    Dk_1[i] = this.Max(matrix[i, k - 1].score - penaltyGapExist, Dk_1[i] - penaltyGapExtend);
                    Dk_1_var[i] = this.Max(matrix[i, k - 1].score - penaltyGapExtend, Dk_1_var[i] - penaltyGapExtend);
                    if (boundaryQ[k - 1] == true)
                    {
                        Dk_1_var[i] = Dk_1_var[i] - boundaryCrossPenalty;

                    }

                    float downScoreInsertion = Math.Max(Dk_1[i], (Dk_1_var[i] + penaltyGapExtend - gp));






                    //Find the highest scroing deletion, testing all matrix entries to the left side
                    Qk_1[k] = this.Max(matrix[i - 1, k].score - penaltyGapExist, Qk_1[k] - penaltyGapExtend);
                    Qk_1_var[k] = this.Max(matrix[i - 1, k].score - penaltyGapExtend, Qk_1_var[k] - penaltyGapExtend);
                    if (boundaryD[i - 1] == true)
                    {
                        Qk_1_var[k] = Qk_1_var[k] - boundaryCrossPenalty;

                    }

                    float rightScoreDeletion = Math.Max(Qk_1[k], (Qk_1_var[k] + penaltyGapExtend - gp));








                    matrix[i, k] = GetMaximumPosition(scoreDiagonal, downScoreInsertion, rightScoreDeletion, 0);

#if DEBUG
                    MatrixPosition mp = GetMaximumPosition(scoreDiagonal, downScoreInsertion, rightScoreDeletion, 0);
#endif



                }
            }
            //Call the traceback method
            Traceback();
        }

        private void ProcessGapExistPenalty()
        {

            float maxStep = Math.Abs(substitutionMatrix.HighestScore) + Math.Abs(substitutionMatrix.LowestScore) - 0.1F;
            float lowestGapExist = penaltyGapExist - maxStep;
            //if the lowestGapExist penalty is smaller than the gapExtend penalty increase it to the gapExtend penalty
            lowestGapExist = (lowestGapExist < penaltyGapExtend ? penaltyGapExtend + 0.1F : lowestGapExist);

            //Initialize the new arrays
            geD = new float[databaseSequence.Length + 1];
            geQ = new float[querySequence.Length + 1];
            boundaryD = new bool[databaseSequence.Length + 1];
            boundaryQ = new bool[querySequence.Length + 1];
            for (int i = 0; i < boundaryD.Length; i++)
            {
                boundaryD[i] = false;
            }
            for (int i = 0; i < boundaryQ.Length; i++)
            {
                boundaryQ[i] = false;
            }




            //The first two characters for every position specific gap penalty matrix are identical to the default gap existence penalty.
            geD[0] = this.penaltyGapExist;
            geQ[0] = this.penaltyGapExist;
            if (databaseSequence.Length >= 1) geD[1] = this.penaltyGapExist;
            if (querySequence.Length >= 1) geQ[1] = this.penaltyGapExist;

            //|0123456789
            //|TTTTTTTTTTTTT
            //99876543211111 geD
            //0123456789

            for (int i = 0; i < databaseSequence.Length - 1; i++)
            {

                if (databaseSequence[i] != databaseSequence[i + 1])
                {
                    geD[i + 2] = penaltyGapExist;
                    boundaryD[i + 1] = true;

                }
                else
                {
                    int k = 1; //k==1 -> unique k==2 AA; k==4 AAAA
                    while (i+k< databaseSequence.Length && databaseSequence[i] == databaseSequence[i + k])
                    {
                        k++;
                    }
                    float step = (penaltyGapExist - lowestGapExist) / (float)(k - 1);

                    //Fill until the new match starts
                    int l = 0;
                    for (; l < k - 1; l++)
                    {
                        geD[i + 2 + l] = geD[i + 1 + l] - step;
#if DEBUG
                        float assign = (geD[i + 1 + l] - step);
                        if (assign <= penaltyGapExtend) throw new Exception("Fucking impossible");
#endif
                    }
                    i = i + l - 1;


                }

            }


            //|0123456789
            //|TTTTTTTTTTTTT
            //99876543211111 geD
            //0123456789
            for (int i = 0; i < querySequence.Length - 1; i++)
            {

                if (querySequence[i] != querySequence[i + 1])
                {
                    geQ[i + 2] = penaltyGapExist;
                    boundaryQ[i + 1] = true;

                }
                else
                {
                    int k = 1; //k==1 -> unique k==2 AA; k==4 AAAA
                    while (i + k < querySequence.Length &&  querySequence[i] == querySequence[i + k])
                    {
                        k++;
                    }
                    float step = (penaltyGapExist - lowestGapExist) / (float)(k - 1);

                    //Fill until the new match starts
                    int l = 0;
                    for (; l < k - 1; l++)
                    {
                        geQ[i + 2 + l] = geQ[i + 1 + l] - step;

#if DEBUG
                        if ((geQ[i + 1 + l] - step) <= penaltyGapExtend) throw new Exception("Fucking impossible");
#endif
                    }

                    i = i + l - 1;


                }

            }





        }


        protected override void Traceback()
        {

            PairwiseAlignmentBuilder builder = new PairwiseAlignmentBuilder();
            //Actual position in the matrix; D..database; Q..query
            int posD = this.databaseSequence.Length;
            int posQ = this.querySequence.Length;
            int lastPosD = posD;
            int lastPosQ = posQ;
            float gp = penaltyGapExist;

            while (matrix[posD, posQ].lastMove != LastMove.None)
            {


                switch (matrix[posD, posQ].lastMove)
                {
                    case LastMove.Diagonal: //while

                        posD--; posQ--;

                        if (posD >= overlap && posQ >= overlap)
                        {
                            lastPosD = posD;
                            lastPosQ = posQ;
                            builder.Append_5_prime(databaseSequence[posD], querySequence[posQ]);
                        }
                        else
                        {
                            goto endOfLoop;
                        }

                        break;
                    case LastMove.Deletion:
                        int startPosD = posD;

                        gp = Math.Min(geD[posD], geQ[posQ]);
                        do
                        {

                            builder.Append_5_prime(databaseSequence[posD - 1], '-');
                            if (gp != penaltyGapExist && boundaryD[posD - 1] == true)
                            {
                                gp = gp + boundaryCrossPenalty;
                                if (gp > penaltyGapExist) gp = penaltyGapExist;
                            }


                            posD--;

                            if (!(posD > overlap)) //eauals posD<=overlap
                            {
                                builder.Remove_5_prime(startPosD - posD);
                                goto endOfLoop;
                            }
#if DEBUG
                            double a = Math.Round(matrix[posD, posQ].score, 4);
                            double b = Math.Round((matrix[startPosD, posQ].score + gp + penaltyGapExtend * (startPosD - posD - 1)), 4);
#endif

                        }
                        while (Math.Abs(Math.Round(matrix[posD, posQ].score, 6) - Math.Round((matrix[startPosD, posQ].score + gp + penaltyGapExtend * (startPosD - posD - 1)), 6)) > 0.0001);

                        break;
                    case LastMove.Insertion:
                        int startPosQ = posQ;
                        gp = Math.Min(geD[posD], geQ[posQ]);
                        do
                        {
                            builder.Append_5_prime('-', querySequence[posQ - 1]);
                            if (gp != penaltyGapExist && boundaryQ[posQ - 1] == true)
                            {
                                gp = gp + boundaryCrossPenalty;
                                if (gp > penaltyGapExist) gp = penaltyGapExist;
                            }


                            posQ--;

                            if (!(posQ > overlap))
                            {
                                builder.Remove_5_prime(startPosQ - posQ);
                                goto endOfLoop;
                            }

#if DEBUG
                            double a = Math.Round(matrix[posD, posQ].score, 4);
                            double b = Math.Round((matrix[posD, startPosQ].score + gp + penaltyGapExtend * (startPosQ - posQ - 1)), 4);
#endif


                        }
                        while (Math.Abs(Math.Round(matrix[posD, posQ].score, 6) - Math.Round((matrix[posD, startPosQ].score + gp + penaltyGapExtend * (startPosQ - posQ - 1)), 6)) > 0.0001);

                        break;
                    default: throw new InvalidOperationException("this should not happen");
                }

            }
        endOfLoop: ;


            if (matrix[databaseSequence.Length, querySequence.Length].score != 0)
            {
                startAlignmentDatabase = -(databaseSequence.Length - lastPosD);
                startAlignmentQuery = -(querySequence.Length - lastPosQ);
                endAlignmentDatabase = -1;
                endAlignmentQuery = -1;
                this.highScore = matrix[databaseSequence.Length, querySequence.Length].score - matrix[lastPosD, lastPosQ].score;
                this.highScoreAtTheBeginning = false;
                alignment = builder.GetAlignment();
            }

            //If the alignment length is zero, reset the values
            if (alignment != null && alignment.Length < 1)
            {
                alignment = null;
                this.highScoreAtTheBeginning = true;
                this.startAlignmentDatabase = 0;
                this.startAlignmentQuery = 0;
                this.endAlignmentDatabase = 0;
                this.endAlignmentQuery = 0;
            }




        }





        public string ToStringVariableGapPenaltyDatabase()
        {

            StringBuilder sb = new StringBuilder("");
            for (int i = 0; i < geD.Length; i++)
            {
                sb.Append((int)geD[i]);
            }
            return sb.ToString();

        }

        public string ToStringVariableGapPenaltyQuery()
        {

            StringBuilder sb = new StringBuilder("");
            for (int i = 0; i < geQ.Length; i++)
            {
                sb.Append((int)geQ[i]);
            }
            return sb.ToString();

        }

        public string ToStringBoundariesDatabase()
        {

            StringBuilder sb = new StringBuilder("");
            for (int i = 0; i < boundaryD.Length; i++)
            {
                sb.Append(boundaryD[i] == true ? "T" : "F");
            }
            return sb.ToString();

        }

        public string ToStringBoundariesQuery()
        {

            StringBuilder sb = new StringBuilder("");
            for (int i = 0; i < boundaryQ.Length; i++)
            {
                sb.Append(boundaryQ[i] == true ? "T" : "F");
            }
            return sb.ToString();

        }






        public bool EndofDynamicExtension
        {
            get
            {
                return this.highScoreAtTheBeginning;
            }
        }


    }



}
