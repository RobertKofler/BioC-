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


namespace Bio.Alignment.Misc
{
    using Bio.Seq;
    using Bio.Alignment;
    using Bio.Seq.IO;

    public abstract class SmithWaterman_Base : IDynamicProgramming
    {
        protected struct Position
        {
            public Position(int database, int query)
            {
                this.database = (ushort)database;
                this.query = (ushort)query;
            }
            public ushort database;
            public ushort query;
        }

        protected struct MatrixPosition
        {
            public MatrixPosition(float score, Position position)
            {
                this.score = score;
                this.position = position;
            }
            public float score;
            public Position position;
        }

        //Core Variables
        protected MatrixPosition[,] matrix = null;
        protected ISequenceContainer databaseSequence;
        protected ISequenceContainer querySequence;
        protected SubstitutionMatrix substitutionMatrix;
        protected PairwiseAlignment alignment = null;


        //Working variables
        protected float highScore = 0;
        protected int endAlignmentDatabase = 0;
        protected int endAlignmentQuery = 0;
        protected int startAlignmentDatabase = 0;
        protected int startAlignmentQuery = 0;


        //Penalties
        protected float penaltyGapExist;
        protected float penaltyGapExtend;

        public SmithWaterman_Base(ISequenceContainer databaseSequence, ISequenceContainer querySequence, SubstitutionMatrix substitutionMatrix)
        {
            if (databaseSequence.Length < 1 || querySequence.Length < 1) throw new InvalidSequenceException("Error in Smith-Waterman dynamic programming. Length of database sequence and query sequence has to be at least 1 bp");

            this.penaltyGapExist = substitutionMatrix.GapExistPenalty;
            this.penaltyGapExtend = substitutionMatrix.GapExtendPenalty;
            this.databaseSequence = databaseSequence;
            this.querySequence = querySequence;
            this.substitutionMatrix = substitutionMatrix;


        }

        protected abstract void InitializeMatrix();




        protected void Traceback()
        {

            PairwiseAlignmentBuilder builder = new PairwiseAlignmentBuilder();

            Position actualPosition = new Position(endAlignmentDatabase, endAlignmentQuery);

            while (matrix[actualPosition.database, actualPosition.query].score != 0)
            {
                Position nextPosition = matrix[actualPosition.database, actualPosition.query].position;
                int deltaQuery = actualPosition.query - nextPosition.query;
                int deltaDatabase = actualPosition.database - nextPosition.database;

                if (deltaQuery == 1 && deltaDatabase == 1)
                {
                    builder.Append_5_prime(databaseSequence[actualPosition.database - 1], querySequence[actualPosition.query - 1]);
                }
                else if (deltaQuery >= 1 && deltaDatabase == 0)
                {
                    //Insertion
                    for (int i = actualPosition.query; i > nextPosition.query; i--)
                    {
                        builder.Append_5_prime('-', querySequence[i - 1]);
                    }
                }
                else if (deltaDatabase >= 1 && deltaQuery == 0)
                {
                    //Deletion
                    for (int i = actualPosition.database; i > nextPosition.database; i--)
                    {
                        builder.Append_5_prime(databaseSequence[i - 1], '-');
                    }
                }
                else throw new Exception("this is soo not good");

                actualPosition = nextPosition;

                /*
                 *                             case LastMove.Insertion:
                    builder.Append_5_prime('-', querySequence[posQ - 1]);
                    posQ--;
                 */
            }



            if (endAlignmentDatabase > 0 && endAlignmentQuery > 0)
            {
                // if (lm != LastMove.Diagonal) throw new Exception("shit this should not happen");
                startAlignmentDatabase = actualPosition.database + 1;
                startAlignmentQuery = actualPosition.query + 1;
                alignment = builder.GetAlignment();

            }



        }



        protected MatrixPosition GetMaximumPosition(float scoreDiagonal, Position positionDiagonal, float scoreInsertion, Position positionInsertion, float scoreDeletion, Position positionDeletion, float scoreNone)
        {
            MatrixPosition position;

            if (scoreDiagonal > scoreNone)
            {
                //exclude scoreNone

                if (scoreDiagonal >= scoreInsertion)
                {
                    //exclude scoreNone & scoreInsertion

                    if (scoreDiagonal >= scoreDeletion)
                    {
                        //exclude scoreNone & scoreInsertion & scoreDeletion => DIAGONAL

                        position = new MatrixPosition(scoreDiagonal, positionDiagonal);
                    }
                    else
                    {
                        //exclude scoreNone & scoreInsertion & scoreDiagonal => DELETION
                        position = new MatrixPosition(scoreDeletion, positionDeletion);

                    }
                }
                else
                {
                    //exclude scoreNone & scoreDiagonal


                    if (scoreInsertion > scoreDeletion)
                    {
                        //exclude scoreNone & scoreDiagonal & scoreDeletion => INSERTION
                        position = new MatrixPosition(scoreInsertion, positionInsertion);
                    }
                    else
                    {
                        //exclude scoreNone &scoreDiagonal & scoreInsertion => DELETION
                        position = new MatrixPosition(scoreDeletion, positionDeletion);

                    }
                }
            }
            else
            {
                //exclude scoreDiagonal
                if (scoreInsertion > scoreNone)
                {
                    //exclude scoreDiagonal & scoreNone

                    if (scoreInsertion > scoreDeletion)
                    {
                        //exclude scoreDiagonal & scoreNone & scoreDeletion => INSERTION
                        position = new MatrixPosition(scoreInsertion, positionInsertion);
                    }
                    else
                    {
                        //exclude scoreDiagonal & scoreNone & scoreInsertion => DELETION
                        position = new MatrixPosition(scoreDeletion, positionDeletion);

                    }
                }
                else
                {
                    //exclude scoreDiagonal & scoreInsertion
                    if (scoreDeletion > scoreNone)
                    {
                        //exclude scoreDiagonal & scoreInsertion & scoreNone => DELETION
                        position = new MatrixPosition(scoreDeletion, positionDeletion);

                    }
                    else
                    {
                        //exclude scoreDiagonal & scoreInsertion & scoreDeletion =>NONE
                        position = new MatrixPosition(scoreNone, new Position(0, 0));
                    }
                }


            }


            return position; //That was annoying
        }


        /// <summary>
        /// Returns the score of the best alignment found for the two specified sequences
        /// </summary>
        public float Score
        {
            get
            {
                return highScore;
            }
        }


        public PairwiseAlignment GetAlignment()
        {

            return alignment;

        }


        public int Start_Database
        {
            get
            {
                return this.startAlignmentDatabase;
            }
        }

        /// <summary>
        /// Start postion of the optimal alignment with respect to the query sequence
        /// </summary>
        public int Start_Query
        {
            get
            {
                return this.startAlignmentQuery;
            }
        }

        /// <summary>
        /// End position of the optimal alignment with respect to the database sequence
        /// </summary>
        public int End_Database
        {
            get
            {
                return this.endAlignmentDatabase;
            }
        }

        /// <summary>
        /// End position of the optimal alignment with respect to the query sequence
        /// </summary>
        public int End_Query
        {
            get
            {
                return this.endAlignmentQuery;
            }
        }


        //Return all values of the matrix in a easy and clearly comprehensible manner
        public override string ToString()
        {


            StringBuilder sb = new StringBuilder("");

            //First write the database sequence
            sb.Append("\t$\t");
            for (int i = 0; i < databaseSequence.Length; i++)
            {
                sb.Append(String.Format("{0}\t", databaseSequence[i]));
            }
            sb.Append('\n');

            for (int k = 0; k < matrix.GetLength(1); k++)
            {
                if (k == 0)
                {
                    sb.Append("$\t");
                }
                else
                {
                    sb.Append(String.Format("{0}\t", querySequence[k - 1]));
                }


                for (int i = 0; i < matrix.GetLength(0); i++)
                {
                    sb.Append(String.Format("{0:f}", matrix[i, k].score));
                    if (i != matrix.GetLength(0) - 1) sb.Append("\t");

                }

                sb.Append("\n");
            }

            return sb.ToString();

        }



        public SubstitutionMatrix SubstitutionMatrix
        {
            get
            {
                return this.substitutionMatrix;
            }


        }

    }



    /// <summary>
    /// Reference implementation of the original Smith Waterman algorithm
    /// This implementation is slow (M^2 * N) and requires much memory
    /// but definitely produces the most reliable results, should only be used for the purpose of testing
    /// uses the equation: Wk= PExi - PExt * (x-1) for calculating the gap penalties
    /// </summary>
    public class SmithWaterman : SmithWaterman_Base
    {


        public SmithWaterman(ISequenceContainer databaseSequence, ISequenceContainer querySequence, SubstitutionMatrix substitutionMatrix)
            : base(databaseSequence, querySequence, substitutionMatrix)
        {


            InitializeMatrix();
        }



        /// <summary>
        /// Initializes the smith waterman alignment matrix, like in the original publication, testing for the highest Wk-insertions
        /// the highest Wk-deletions and than max(s(i,k),wk-del,wk-ins,0);
        /// </summary>
        protected override void InitializeMatrix()
        {
            matrix = new MatrixPosition[databaseSequence.Length + 1, querySequence.Length + 1];

            //      @   A   C   A   G   Database
            //  @   0   0   0   0   0
            //  Q   0
            //  U   0
            //  E   0
            //  R   0
            //  Y   0
            //  
            //DATABASE  AATTCC-GGG a '-' in the database is a insertion, a insertion is a move down
            //QUERY     AA-TCCGGGG a '-' in the query is a deletion, a deletion is a move to the right


            //First set all values in the horizontal to zero
            for (int i = 0; i < matrix.GetLength(0); i++)
            {
                matrix[i, 0] = new MatrixPosition(0, new Position(0, 0));
            }

            //Second set all values in the vertical to zero
            for (int i = 0; i < matrix.GetLength(1); i++)
            {
                matrix[0, i] = new MatrixPosition(0, new Position(0, 0));
            }

            //Go down use dimension of the query
            for (int k = 1; k < matrix.GetLength(1); k++)
            {
                for (int i = 1; i < matrix.GetLength(0); i++)
                {
                    //i=database sequence in the horizontal
                    //k=query sequence in the vertical

                    Position bestDeletionPosition = new Position(0, 0);
                    Position bestInsertionPosition = new Position(0, 0);

                    //the database sequence is in the horizontal, the query in the vertical axis of the matrix
                    //Diagonal score is the previous score and in addition the similarityValue;
                    float scoreDiagonal = matrix[i - 1, k - 1].score + substitutionMatrix.GetSimilarityValue(databaseSequence[i - 1], querySequence[k - 1]);

                    //Find the highest scoring insertion, testing all matrix to the upper side;
                    float downScoreInsertion = 0;
                    for (int down = k - 1; down >= 0; down--)
                    {
                        if (matrix[i, down].score > scoreDiagonal + penaltyGapExist)
                        {
                            float attempt = matrix[i, down].score - penaltyGapExist - penaltyGapExtend * (k - down - 1);

                            //Found a better insertion
                            if (attempt > downScoreInsertion)
                            {
                                downScoreInsertion = attempt;
                                bestInsertionPosition = new Position(i, down);
                            }
                        }

                    }

                    //Find the highest scroing deletion, testing all matrix entries to the left side
                    float rightScoreDeletion = 0;
                    for (int left = i - 1; left >= 0; left--)
                    {
                        if (matrix[left, k].score > scoreDiagonal + penaltyGapExist)
                        {
                            float attempt = matrix[left, k].score - penaltyGapExist - penaltyGapExtend * (i - left - 1);
                            //Found a better deletion
                            if (attempt > rightScoreDeletion)
                            {
                                rightScoreDeletion = attempt;
                                bestDeletionPosition = new Position(left, k);
                            }
                        }

                    }

                    matrix[i, k] = this.GetMaximumPosition(scoreDiagonal, new Position(i - 1, k - 1), downScoreInsertion, bestInsertionPosition, rightScoreDeletion, bestDeletionPosition, 0);

                    //Updating the highest scoring matrix entry
                    if (matrix[i, k].score > highScore)
                    {
                        //new highscore
                        highScore = matrix[i, k].score;
                        endAlignmentDatabase = i;
                        endAlignmentQuery = k;
                    }



                }
            }
            Traceback();
        }


    }


    /// <summary>
    /// Reference implementation of the original Smith Waterman algorithm
    /// This implementation is slow (M^2 * N) and requires much memory
    /// but definitely produces the most reliable results, should only be used for the purpose of testing
    /// special adaption for 454 sequencing, uses a decreasing gap exist penalty and boundaries
    /// </summary>
    public class SmithWaterman_454P : SmithWaterman_Base
    {


        //454 Adaption
        private float[] geD; //gap exist database
        private bool[] boundaryD;
        private bool[] boundaryQ;
        private float[] geQ; //gap exist query
        private float gapLength = 6.0F;


        public SmithWaterman_454P(ISequenceContainer databaseSequence, ISequenceContainer querySequence, SubstitutionMatrix substitutionMatrix)
            : base(databaseSequence, querySequence, substitutionMatrix)
        {

            InitializeMatrix();
        }



        /// <summary>
        /// Initializes the smith waterman alignment matrix, like in the original publication, testing for the highest Wk-insertions
        /// the highest Wk-deletions and than max(s(i,k),wk-del,wk-ins,0);
        /// </summary>
        protected override void InitializeMatrix()
        {
            matrix = new MatrixPosition[databaseSequence.Length + 1, querySequence.Length + 1];
            ProcessGapExistPenalty();

            //      @   A   C   A   G   Database
            //  @   0   0   0   0   0
            //  Q   0
            //  U   0
            //  E   0
            //  R   0
            //  Y   0
            //  
            //DATABASE  AATTCC-GGG a '-' in the database is a insertion, a insertion is a move down
            //QUERY     AA-TCCGGGG a '-' in the query is a deletion, a deletion is a move to the right


            //First set all values in the horizontal to zero
            for (int i = 0; i < matrix.GetLength(0); i++)
            {
                matrix[i, 0] = new MatrixPosition(0, new Position(0, 0));
            }

            //Second set all values in the vertical to zero
            for (int i = 0; i < matrix.GetLength(1); i++)
            {
                matrix[0, i] = new MatrixPosition(0, new Position(0, 0));
            }

            //Go down use dimension of the query
            for (int k = 1; k < matrix.GetLength(1); k++)
            {
                for (int i = 1; i < matrix.GetLength(0); i++)
                {
                    //i=database sequence in the horizontal
                    //k=query sequence in the vertical
                    float gp = Math.Min(geD[i], geQ[k]);

                    Position bestDeletionPosition = new Position(0, 0);
                    Position bestInsertionPosition = new Position(0, 0);

                    //the database sequence is in the horizontal, the query in the vertical axis of the matrix
                    //Diagonal score is the previous score and in addition the similarityValue;
                    float scoreDiagonal = matrix[i - 1, k - 1].score + substitutionMatrix.GetSimilarityValue(databaseSequence[i - 1], querySequence[k - 1]);

                    //Find the highest scoring insertion, testing all matrix to the upper side;
                    float downScoreInsertion = 0;
                    bool crossBoundary = false;
                    for (int down = k - 1; down >= 0; down--)
                    {
                        if (matrix[i, down].score > scoreDiagonal)
                        {

                            //Test if a boundary is crossed
                            if (boundaryQ[down] == true) crossBoundary = true;
                            if (crossBoundary) gp = penaltyGapExist;


                            float attempt = matrix[i, down].score - gp - penaltyGapExtend * (k - down - 1);

                            //Found a better insertion
                            if (attempt > downScoreInsertion)
                            {
                                downScoreInsertion = attempt;
                                bestInsertionPosition = new Position(i, down);
                            }
                        }
                        else if (matrix[i, down].score == 0) break;

                    }

                    //Find the highest scroing deletion, testing all matrix entries to the left side

                    gp = Math.Min(geD[i], geQ[k]);
                    float rightScoreDeletion = 0;
                    crossBoundary = false;
                    for (int left = i - 1; left >= 0; left--)
                    {
                        if (matrix[left, k].score > scoreDiagonal)
                        {

                            //Test if a boundary is crossed
                            if (boundaryD[left] == true) crossBoundary = true;
                            if (crossBoundary) gp = penaltyGapExist;


                            float attempt = matrix[left, k].score - gp - penaltyGapExtend * (i - left - 1);
                            //Found a better deletion
                            if (attempt > rightScoreDeletion)
                            {
                                rightScoreDeletion = attempt;
                                bestDeletionPosition = new Position(left, k);
                            }
                        }
                        else if (matrix[left, k].score == 0) break;

                    }

                    matrix[i, k] = this.GetMaximumPosition(scoreDiagonal, new Position(i - 1, k - 1), downScoreInsertion, bestInsertionPosition, rightScoreDeletion, bestDeletionPosition, 0);

                    //Updating the highest scoring matrix entry
                    if (matrix[i, k].score > highScore)
                    {
                        //new highscore
                        highScore = matrix[i, k].score;
                        endAlignmentDatabase = i;
                        endAlignmentQuery = k;
                    }



                }
            }
            Traceback();
        }

        private void ProcessGapExistPenalty()
        {
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
            float defaultStep = (penaltyGapExist - penaltyGapExtend) / gapLength;

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
                    else step = (penaltyGapExist - penaltyGapExtend) / (float)k;

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
                    //if (k > 1) boundaryD[i + 1 + l] = true;
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
                    else step = (penaltyGapExist - penaltyGapExtend) / (float)k;

                    //Fill until the new match starts
                    int l = 0;
                    for (; l < k - 1; l++)
                    {
                        geQ[i + 2 + l] = geQ[i + 1 + l] - step;

#if DEBUG
                        if ((geQ[i + 1 + l] - step) <= penaltyGapExtend) throw new Exception("Fucking impossible");
#endif
                    }
                    //if (k > 1) boundaryQ[i + 1 + l] = true;

                    i = i + l - 1;


                }

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
    





