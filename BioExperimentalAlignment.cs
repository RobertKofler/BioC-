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


namespace Bio.ExperimentalAlignment
{
    using Bio.Seq;
    using Bio.Alignment;
    using Bio.Seq.IO;

    /// <summary>
    /// Implementation of the Smith Waterman algorithm, modified according to Gotoh.
    /// Specialy adapted for the 454 sequencing, uses boundarie cross penaltys for gaps extending beyond poly-N tracts
    /// The direction of the 454 read is important, the query sequence has to have the 5'->3' direction 
    /// P..Plus/Plus
    /// </summary>
    public class NeedlemanWunschGotoh_454P : SmithWatermanGotoh
    {
        private float[] geD; //gap exist database
        private float[] geQ; //gap exist query
        private bool[] boundaryD;
        private bool[] boundaryQ;

        private float[] Dk_1_var;
        private bool[] dk_1_var_boundary;
        private bool[] qk_1_var_boundary;
        private float[] Qk_1_var;
        private float gapLength = 3.0F;
        private float boundaryCrossPenalty;

        protected NeedlemanWunschGotoh_454P()
        {
        }

        public NeedlemanWunschGotoh_454P(ISequenceContainer databaseSequence, ISequenceContainer querySequence, SubstitutionMatrix substitutionMatrix)
        {
            if (databaseSequence.Length < 1 || querySequence.Length < 1) throw new InvalidSequenceException("Error in Smith-Waterman-Gotoh dynamic programming. Length of database sequence or query sequence has to be at least 1 bp");


            this.penaltyGapExist = substitutionMatrix.GapExistPenalty;
            this.penaltyGapExtend = substitutionMatrix.GapExtendPenalty;
            this.databaseSequence = databaseSequence;
            this.querySequence = querySequence;
            this.substitutionMatrix = substitutionMatrix;

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
            Qk_1[0] = Qk_1_var[0] = 0;
            Dk_1[0] = Dk_1_var[0] = 0;

            ProcessGapExistPenalty(); //Create position specific gap existence penalty matrix

            //First set all values in the horizontal (database) to zero
            int boundary = 0;
            for (int i = 1; i < matrix.GetLength(0); i++)
            {
                if (boundaryD[i - 1] == true) boundary++;
                float gap_fix = 0.0F - penaltyGapExist - (penaltyGapExtend * (i - 1)); //merely the gap extending penalty
                float gp = Math.Min(geD[i], geQ[0]);
                float gap_var = 0.0F - boundaryCrossPenalty * boundary - (penaltyGapExtend * (i - 1)) - gp;
                matrix[i, 0] = new MatrixPosition(Math.Max(gap_fix, gap_var), LastMove.Deletion);
                Dk_1[i] = matrix[i, 0].score - penaltyGapExist;
                Dk_1_var[i] = matrix[i, 0].score - penaltyGapExtend;
                if (boundaryQ[0] == true) Dk_1_var[i] -= boundaryCrossPenalty;


            }
            //Second set all values in the vertical (query) to zero
            boundary = 0;
            for (int i = 1; i < matrix.GetLength(1); i++)
            {
                if (boundaryQ[i - 1] == true) boundary++;
                float gap_fix = 0.0F - penaltyGapExist - (penaltyGapExtend * (i - 1)); //merely the gap extending penalty
                float gp = Math.Min(geD[0], geQ[i]);
                float gap_var = 0.0F - boundaryCrossPenalty * boundary - (penaltyGapExtend * (i - 1)) - gp;
                matrix[0, i] = new MatrixPosition(Math.Max(gap_fix, gap_var), LastMove.Insertion);

                Qk_1[i] = matrix[0, i].score - penaltyGapExist;
                Qk_1_var[i] = matrix[0, i].score - penaltyGapExtend;
                if (boundaryD[0] == true) Qk_1_var[i] -= boundaryCrossPenalty;

            }



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


                }
            }

            //Call the traceback method
            Traceback();
        }

        private void ProcessGapExistPenalty()
        {

            //Roberts Law
            float lowestGapExist = 2 * penaltyGapExtend; // +0.1F; //substitutionMatrix.LowestScore should be negative
            // float lowestGapExist = penaltyGapExist - substitutionMatrix.HighestScore + substitutionMatrix.LowestScore + 0.1F; //substitutionMatrix.LowestScore should be negative


            //if the lowestGapExist penalty is smaller than the gapExtend penalty increase it to the gapExtend penalty
            lowestGapExist = (lowestGapExist < penaltyGapExtend ? penaltyGapExtend + 0.1F : lowestGapExist);


            //CRossboundary < Hit +mmp
            //Roberts law for boundaries Lp > Lp+ CrossBoundary - hit - mismatchpenalty -> crossboundary < hit+ mismatchpenalty
            this.boundaryCrossPenalty = 2.0F * penaltyGapExtend;
            if (boundaryCrossPenalty >= substitutionMatrix.HighestScore + Math.Abs(substitutionMatrix.LowestScore)) boundaryCrossPenalty = substitutionMatrix.HighestScore + Math.Abs(substitutionMatrix.LowestScore) - 0.1F;


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


            //The decreas of the gapExist penalty with each step
            float defaultStep = (penaltyGapExist - lowestGapExist) / (gapLength - 1);

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
                    while (databaseSequence[i] == databaseSequence[i + k])
                    {
                        k++;
                    }
                    float step;

                    //Set the steps
                    if (k < gapLength) step = defaultStep;
                    else step = (penaltyGapExist - lowestGapExist) / (float)(k - 1);

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
                    while (querySequence[i] == querySequence[i + k])
                    {
                        k++;
                    }
                    float step;

                    //Set the steps
                    if (k < gapLength) step = defaultStep;
                    else step = (penaltyGapExist - lowestGapExist) / (float)(k - 1);

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
            int posD = databaseSequence.Length;
            int posQ = querySequence.Length;
            float gp = penaltyGapExist;

            while (posD != 0 || posQ != 0)
            {

                bool crossBoundary = false;

                switch (matrix[posD, posQ].lastMove)
                {
                    case LastMove.Diagonal: //while
                        builder.Append_5_prime(databaseSequence[posD - 1], querySequence[posQ - 1]);
                        posD--; posQ--;

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
                            // if (posD == 0 && posQ == 0) break;
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
                            //if (posD == 0 && posQ == 0) break;

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

            if (databaseSequence.Length > 0 && querySequence.Length > 0)
            {
                startAlignmentDatabase = 1;
                startAlignmentQuery = 1;
                endAlignmentDatabase = databaseSequence.Length;
                endAlignmentQuery = querySequence.Length;
                alignment = builder.GetAlignment();
                this.highScore = matrix[endAlignmentDatabase, endAlignmentQuery].score;
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



    }

    /// <summary>
    /// Implementation of the Smith Waterman algorithm, modified according to Gotoh.
    /// Specialy adapted for the 454 sequencing, uses boundarie cross penaltys for gaps extending beyond poly-N tracts
    /// The direction of the 454 read is important, the query sequence has to have the 5'->3' direction 
    /// P..Plus/Plus
    /// </summary>
    public class NeedlemanWunschGotoh_454M : SmithWatermanGotoh
    {
        private float[] geD; //gap exist database
        private float[] geQ; //gap exist query
        private bool[] boundaryD;
        private bool[] boundaryQ;

        private float[] Dk_1_var;
        private bool[] dk_1_var_boundary;
        private bool[] qk_1_var_boundary;
        private float[] Qk_1_var;
        private float gapLength = 3.0F;
        private float boundaryCrossPenalty;



        protected NeedlemanWunschGotoh_454M()
        {
        }

        public NeedlemanWunschGotoh_454M(ISequenceContainer databaseSequence, ISequenceContainer querySequence, SubstitutionMatrix substitutionMatrix)
        {
            if (databaseSequence.Length < 1 || querySequence.Length < 1) throw new InvalidSequenceException("Error in Smith-Waterman-Gotoh dynamic programming. Length of database sequence or query sequence has to be at least 1 bp");


            this.penaltyGapExist = substitutionMatrix.GapExistPenalty;
            this.penaltyGapExtend = substitutionMatrix.GapExtendPenalty;
            this.databaseSequence = databaseSequence;
            this.querySequence = querySequence;
            this.substitutionMatrix = substitutionMatrix;



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

            this.Qk_1[querySequence.Length] = 0;
            this.Dk_1[databaseSequence.Length] = 0;
            Qk_1_var[querySequence.Length] = 0;
            Dk_1_var[databaseSequence.Length] = 0;

            ProcessGapExistPenalty(); //Create position specific gap existence penalty matrix


            //First set all values in the horizontal (database) to zero
            int boundary = 0;
            for (int i = databaseSequence.Length - 1; i >= 0; i--)
            {
                int dist = databaseSequence.Length - i;
                if (boundaryD[i + 1] == true) boundary++;
                float gap_fix = 0.0F - penaltyGapExist - (penaltyGapExtend * (dist - 1)); //merely the gap extending penalty
                float gp = Math.Min(geD[i], geQ[querySequence.Length]);
                float gap_var = 0.0F - boundaryCrossPenalty * boundary - (penaltyGapExtend * (dist - 1)) - gp;
                matrix[i, querySequence.Length] = new MatrixPosition(Math.Max(gap_fix, gap_var), LastMove.Deletion);
                Dk_1[i] = matrix[i, querySequence.Length].score - penaltyGapExist;
                Dk_1_var[i] = matrix[i, querySequence.Length].score - penaltyGapExtend;
                if (boundaryQ[querySequence.Length] == true) Dk_1_var[i] -= boundaryCrossPenalty;
            }

            //Second set all values in the vertical (query) to zero
            boundary = 0;
            for (int i = querySequence.Length - 1; i >= 0; i--)
            {
                int dist = querySequence.Length - i;
                if (boundaryQ[i + 1] == true) boundary++;
                float gap_fix = 0.0F - penaltyGapExist - (penaltyGapExtend * (dist - 1)); //merely the gap extending penalty
                float gp = Math.Min(geD[databaseSequence.Length], geQ[i]);
                float gap_var = 0.0F - boundaryCrossPenalty * boundary - (penaltyGapExtend * (dist - 1)) - gp;
                matrix[databaseSequence.Length, i] = new MatrixPosition(Math.Max(gap_fix, gap_var), LastMove.Deletion);
                Qk_1[i] = matrix[databaseSequence.Length, i].score - penaltyGapExist;
                Qk_1_var[i] = matrix[databaseSequence.Length, i].score - penaltyGapExtend;
                if (boundaryD[databaseSequence.Length] == true) Qk_1_var[i] -= boundaryCrossPenalty;
            }


            //Go down use dimension of the query
            for (int k = querySequence.Length - 1; k >= 0; k--)
            {
                for (int i = databaseSequence.Length - 1; i >= 0; i--)
                {
                    float gp = Math.Min(geD[i], geQ[k]);

                    float scoreDiagonal = matrix[i + 1, k + 1].score + substitutionMatrix.GetSimilarityValue(databaseSequence[i], querySequence[k]);


                    //Find the highest scoring insertion, testing all matrix to the upper side;
                    Dk_1[i] = this.Max(matrix[i, k + 1].score - penaltyGapExist, Dk_1[i] - penaltyGapExtend);
                    Dk_1_var[i] = this.Max(matrix[i, k + 1].score - penaltyGapExtend, Dk_1_var[i] - penaltyGapExtend);
                    if (boundaryQ[k + 1] == true)
                    {
                        Dk_1_var[i] = Dk_1_var[i] - boundaryCrossPenalty;
                    }
                    float downScoreInsertion = Math.Max(Dk_1[i], (Dk_1_var[i] + penaltyGapExtend - gp));

                    //Find the highest scroing deletion, testing all matrix entries to the left side
                    Qk_1[k] = this.Max(matrix[i + 1, k].score - penaltyGapExist, Qk_1[k] - penaltyGapExtend);
                    Qk_1_var[k] = this.Max(matrix[i + 1, k].score - penaltyGapExtend, Qk_1_var[k] - penaltyGapExtend);
                    if (boundaryD[i + 1] == true)
                    {
                        Qk_1_var[k] = Qk_1_var[k] - boundaryCrossPenalty;

                    }
                    float rightScoreDeletion = Math.Max(Qk_1[k], (Qk_1_var[k] + penaltyGapExtend - gp));
                    matrix[i, k] = GetMaximumPosition(scoreDiagonal, downScoreInsertion, rightScoreDeletion);

#if DEBUG
                    MatrixPosition mp = GetMaximumPosition(scoreDiagonal, downScoreInsertion, rightScoreDeletion);
#endif



                }
            }
            //Call the traceback method
            Traceback();
        }

        private void ProcessGapExistPenalty()
        {

            //Roberts Law
            float lowestGapExist = 2 * penaltyGapExtend; // +0.1F; //substitutionMatrix.LowestScore should be negative
            // float lowestGapExist = penaltyGapExist - substitutionMatrix.HighestScore + substitutionMatrix.LowestScore + 0.1F; //substitutionMatrix.LowestScore should be negative


            //if the lowestGapExist penalty is smaller than the gapExtend penalty increase it to the gapExtend penalty
            lowestGapExist = (lowestGapExist < penaltyGapExtend ? penaltyGapExtend + 0.1F : lowestGapExist);


            //CRossboundary < Hit +mmp
            //Roberts law for boundaries Lp > Lp+ CrossBoundary - hit - mismatchpenalty -> crossboundary < hit+ mismatchpenalty
            this.boundaryCrossPenalty = 2.0F * penaltyGapExtend;
            if (boundaryCrossPenalty >= substitutionMatrix.HighestScore + Math.Abs(substitutionMatrix.LowestScore)) boundaryCrossPenalty = substitutionMatrix.HighestScore + Math.Abs(substitutionMatrix.LowestScore) - 0.1F;


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


            //The decreas of the gapExist penalty with each step
            float defaultStep = (penaltyGapExist - lowestGapExist) / (gapLength - 1);

            //The first two characters for every position specific gap penalty matrix are identical to the default gap existence penalty.
            geD[databaseSequence.Length] = this.penaltyGapExist;
            geQ[querySequence.Length] = this.penaltyGapExist;
            if (databaseSequence.Length >= 1) geD[databaseSequence.Length - 1] = this.penaltyGapExist;
            if (querySequence.Length >= 1) geQ[querySequence.Length - 1] = this.penaltyGapExist;
            //01234556789
            // CGTCGTCGT|

            for (int i = databaseSequence.Length - 1; i > 0; i--)
            {

                if (databaseSequence[i] != databaseSequence[i - 1])
                {
                    geD[i - 1] = penaltyGapExist;
                    boundaryD[i] = true;

                }
                else
                {
                    int k = 1; //k==1 -> unique k==2 AA; k==4 AAAA
                    while (i - k >= 0 && databaseSequence[i] == databaseSequence[i - k])
                    {
                        k++;
                    }
                    float step;

                    //Set the steps
                    if (k < gapLength) step = defaultStep;
                    else step = (penaltyGapExist - lowestGapExist) / (float)(k - 1);

                    //Fill until the new match starts
                    int l = 0;
                    for (; l < k - 1; l++)
                    {
                        geD[i - 1 - l] = geD[i - l] - step;
#if DEBUG
                        float assign = (geD[i - l] - step);
                        if (assign <= penaltyGapExtend) throw new Exception("Fucking impossible");
#endif
                    }
                    i = i - l + 1;// i-l+1?


                }

            }


            //|0123456789
            //|TTTTTTTTTTTTT
            //99876543211111 geD
            //0123456789
            for (int i = querySequence.Length - 1; i > 0; i--)
            {

                if (querySequence[i] != querySequence[i - 1])
                {
                    geQ[i - 1] = penaltyGapExist;
                    boundaryQ[i] = true;

                }
                else
                {
                    int k = 1; //k==1 -> unique k==2 AA; k==4 AAAA
                    while (i - k >= 0 && querySequence[i] == querySequence[i - k])
                    {
                        k++;
                    }
                    float step;

                    //Set the steps
                    if (k < gapLength) step = defaultStep;
                    else step = (penaltyGapExist - lowestGapExist) / (float)(k - 1);

                    //Fill until the new match starts
                    int l = 0;
                    for (; l < k - 1; l++)
                    {
                        geQ[i - 1 - l] = geQ[i - l] - step;

#if DEBUG
                        if ((geQ[i - l] - step) <= penaltyGapExtend) throw new Exception("Fucking impossible");
#endif
                    }

                    i = i - l + 1;


                }

            }





        }


        protected override void Traceback()
        {

            PairwiseAlignmentBuilder builder = new PairwiseAlignmentBuilder();
            //Actual position in the matrix; D..database; Q..query
            int posD = 0;
            int posQ = 0;
            float gp = penaltyGapExist;

            while (posD != databaseSequence.Length || posQ != querySequence.Length)
            {

                // bool crossBoundary = false;

                switch (matrix[posD, posQ].lastMove)
                {
                    case LastMove.Diagonal: //while
                        builder.Append_3_prime(databaseSequence[posD], querySequence[posQ]);
                        posD++; posQ++;

                        break;
                    case LastMove.Deletion:
                        int startPosD = posD;

                        gp = Math.Min(geD[posD], geQ[posQ]);
                        do
                        {

                            builder.Append_3_prime(databaseSequence[posD], '-');
                            if (gp != penaltyGapExist && boundaryD[posD + 1] == true)
                            {
                                gp = gp + boundaryCrossPenalty;
                                if (gp > penaltyGapExist) gp = penaltyGapExist;
                            }


                            posD++;
#if DEBUG
                            double a = Math.Round(matrix[posD, posQ].score, 4);
                            double b = Math.Round((matrix[startPosD, posQ].score + gp + penaltyGapExtend * (posD - startPosD - 1)), 4);
#endif

                        }
                        while (Math.Abs(Math.Round(matrix[posD, posQ].score, 6) - Math.Round((matrix[startPosD, posQ].score + gp + penaltyGapExtend * (posD - startPosD - 1)), 6)) > 0.0001);

                        break;
                    case LastMove.Insertion:
                        int startPosQ = posQ;
                        gp = Math.Min(geD[posD], geQ[posQ]);
                        do
                        {
                            builder.Append_3_prime('-', querySequence[posQ]);
                            if (gp != penaltyGapExist && boundaryQ[posQ + 1] == true)
                            {
                                gp = gp + boundaryCrossPenalty;
                                if (gp > penaltyGapExist) gp = penaltyGapExist;
                            }


                            posQ++;

#if DEBUG
                            double a = Math.Round(matrix[posD, posQ].score, 4);
                            double b = Math.Round((matrix[posD, startPosQ].score + gp + penaltyGapExtend * (posQ - startPosQ - 1)), 4);
#endif


                        }
                        while (Math.Abs(Math.Round(matrix[posD, posQ].score, 6) - Math.Round((matrix[posD, startPosQ].score + gp + penaltyGapExtend * (posQ - startPosQ - 1)), 6)) > 0.0001);

                        break;
                    default: throw new InvalidOperationException("this should not happen");
                }


            }
            if (databaseSequence.Length > 0 && querySequence.Length > 0)
            {
                endAlignmentDatabase = databaseSequence.Length;
                endAlignmentQuery = querySequence.Length;
                alignment = builder.GetAlignment();
                startAlignmentQuery = 1;
                startAlignmentDatabase = 1;
                this.highScore = matrix[0, 0].score;
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



    }



}