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

namespace Bio.Alignment
{
    using Bio.Seq;

    /// <summary>
    /// Represents a nucleotide or amino acid substitution matrix.
    /// Substitution matrices are used in Smith-Waterman, Needleman-Wunsch algorithm, phylogenetic studies, BLAST etc.
    /// Can be constructed for amino acids, nucleotides etc. It is possible to define similarity and distance values.
    /// This implementation also contains the gapExist and gapExtend penalties.
    /// </summary>
    public class SubstitutionMatrix
    {
        private Dictionary<string, float> matrix;

        private bool errorWhenNotInMatrix = true;
        private float? gapExistPenalty = null;
        private float? gapExtendPenalty = null;
        private float? highestScore = null;
        private float? lowestScore = null;





        public SubstitutionMatrix(Dictionary<string, float> matrix)
        {
            this.matrix = matrix;
        }

        /// <summary>
        /// Instantiates a representation of a substitution matrix with its matrix;
        /// Returns a default score if no matrix entry is available for the specified key pair (database-char vs query-char)
        /// For instance, this might be necessary for characters like YRWSMK etc in nucleotide sequences, since many matrixes will not explicitely assign values to this rare character combinations
        /// </summary>
        /// <param name="matrix">the substitution matrix in form of a hash table</param>
        ///<param name="errorIfNotInMatrix">should the substitution matrix throw an error if the matrix does not contain a given pair of characters, otherwise the lowest score will be returned</param>
        public SubstitutionMatrix(Dictionary<string, float> matrix, bool errorIfNotInMatrix)
            : this(matrix)
        {

            errorWhenNotInMatrix = errorIfNotInMatrix;
        }

        /// <summary>
        /// Instantiates a representation of a substitution matrix with its matrix;
        /// Returns a default score if no matrix entry is available for the specified key pair (database-char vs query-char)
        /// For instance, this might be necessary for characters like YRWSMK etc in nucleotide sequences, since many matrixes will not explicitely assign values to this rare character combinations
        /// Specify positive values only
        /// </summary>
        /// <param name="matrix">the substitution matrix in form of a hash table</param>
        ///<param name="errorIfNotInMatrix">should the substitution matrix throw an error if the matrix does not contain a given pair of characters, otherwise the lowest score will be returned</param>
        public SubstitutionMatrix(Dictionary<string, float> matrix, bool errorIfNotInMatrix, float gapExistPenalty, float gapExtendPenalty)
            : this(matrix, errorIfNotInMatrix)
        {
            this.gapExistPenalty = gapExistPenalty;
            this.gapExtendPenalty = gapExtendPenalty;
        }

        /// <summary>
        /// Instantiates a representation of a substitution matrix with its matrix;
        /// Returns a default score if no matrix entry is available for the specified key pair (database-char vs query-char)
        /// For instance, this might be necessary for characters like YRWSMK etc in nucleotide sequences, since many matrixes will not explicitely assign values to this rare character combinations
        /// Specify positive values only
        /// </summary>
        /// <param name="matrix">the substitution matrix in form of a hash table</param>
        ///<param name="errorIfNotInMatrix">should the substitution matrix throw an error if the matrix does not contain a given pair of characters, otherwise the lowest score will be returned</param>
        public SubstitutionMatrix(Dictionary<string, float> matrix, bool errorIfNotInMatrix, float gapExistPenalty, float gapExtendPenalty, float highestScore, float lowestScore)
            : this(matrix, errorIfNotInMatrix, gapExistPenalty, gapExtendPenalty)
        {

            this.lowestScore = lowestScore;
            this.highestScore = highestScore;
        }


        /// <summary>
        /// Returns the value (either positive or negative) reflecting the similarity of the two specified characters
        /// </summary>
        public float this[char databaseChar, char queryChar]
        {
            get
            {
                return GetSimilarityValue(databaseChar, queryChar);
            }
        }
        /// <summary>
        /// Returns the value (either positiv or negativ) reflecting the similarity of the two specified characters
        /// </summary>
        public float GetSimilarityValue(char databaseChar, char queryChar)
        {
            string key = databaseChar.ToString() + queryChar.ToString();
            if (matrix.ContainsKey(key))
            {
                return matrix[key];
            }
            else if (!errorWhenNotInMatrix)
            {
                return this.lowestScore.Value;
            }
            else throw new InvalidOperationException("Substitution matrix does not contain the specified characters; Try to either get new substitution matrix or specify default penalty for characters not present in the matrix");
        }

        /// <summary>
        /// Returns the gap exist penalty
        /// </summary>
        public float GapExistPenalty
        {
            get
            {
                return gapExistPenalty.Value;

            }

        }


        /// <summary>
        /// Returns the gap extending penalty
        /// </summary>
        public float GapExtendPenalty
        {
            get
            {
                return gapExtendPenalty.Value;
            }

        }

        /// <summary>
        /// Returns the highest score of the substitution matrix
        /// </summary>
        public float HighestScore
        {
            get
            {
                return this.highestScore.Value;
            }

        }

        /// <summary>
        /// Returns the lowest score of the substitution matrix, should be negative if smith-waterman algorithm will be used
        /// </summary>
        public float LowestScore
        {
            get
            {
                return this.lowestScore.Value;
            }
        }

        /// <summary>
        /// Recalculate the score for a given alignment using the specified hit-scores, mismatch penalties gapexist penalties and gapextend penalties
        /// Careful when recalculating the values for alignments containing non afine gap penalties or alignments which use variing gapExist penalties
        /// </summary>
        /// <param name="alignment"></param>
        /// <returns></returns>
        public float GetScoreForAlignment(PairwiseAlignment alignment)
        {
            float score = 0;
            ISequenceContainer data = alignment.DatabaseSequence;
            ISequenceContainer query = alignment.QuerySequence;

            for (int i = 0; i < alignment.Length; i++)
            {
                int gapCount = 0;
                while (data[i + gapCount] == '-') gapCount++;

                if (gapCount > 0)
                {
                    score -= (this.GapExistPenalty + ((gapCount - 1) * this.GapExtendPenalty));
                    i = i + gapCount - 1;
                    continue;
                }

                gapCount = 0;
                while (query[i + gapCount] == '-') gapCount++;

                if (gapCount > 0)
                {
                    score -= (this.GapExistPenalty + ((gapCount - 1) * this.GapExtendPenalty));
                    i = i + gapCount - 1;
                    continue;
                }

                score += this.GetSimilarityValue(data[i], query[i]);


            }

            return score;
        }

    }


    /// <summary>
    /// Static factory for retrieval of substitution matrices, whether amino acids or nucleotide substitutions 
    /// </summary>
    public static class SubstitutionMatrixFactory
    {
        private static float scoreNX = 0;

        /// <summary>
        /// Get the default substitution matrix for nucleotide sequences, a hit-score and a mismatch-penalty have to be specified
        /// (Recommended scoreHit=1, mismatchPenaltyBaseSubstituion=-1;
        /// </summary>
        /// <param name="scoreHit">specifiy a positive number, which will represent the log ods score for a nucleotide sequence hit, identity</param>
        /// <param name="mismatchPenaltyBaseSubstitution">specify a positve number which will represent the mismatch penalty when two characters are not identical, the program uses the negative value of the specified number</param>
        /// <param name="gapExistPenalty">specify a positive number which represents the gap exist penalty</param>
        /// <param name="gapExtendPenalty">specify a positive number which represents the gap extend penalty</param>
        /// <returns>a substitution matrix having the specified characteristics, can be used in Smith-Waterman or Needleman-Wunsch implementations</returns>
        public static SubstitutionMatrix GetNucleotideSequenceSubstitutionMatrix(float scoreHit, float mismatchPenaltyBaseSubstitution, float gapExistPenalty, float gapExtendPenalty)
        {
            if (scoreHit <= 0 || mismatchPenaltyBaseSubstitution <= 0) throw new ArgumentOutOfRangeException("Please specify meangingfull hit-scores and mismatch penalties for the substitution matrix");
            Dictionary<string, float> matrix = new Dictionary<string, float>();
            matrix.Add("AA", scoreHit);
            matrix.Add("AT", -mismatchPenaltyBaseSubstitution);
            matrix.Add("AC", -mismatchPenaltyBaseSubstitution);
            matrix.Add("AG", -mismatchPenaltyBaseSubstitution);
            matrix.Add("AN", scoreNX);

            matrix.Add("TA", -mismatchPenaltyBaseSubstitution);
            matrix.Add("TT", scoreHit);
            matrix.Add("TC", -mismatchPenaltyBaseSubstitution);
            matrix.Add("TG", -mismatchPenaltyBaseSubstitution);
            matrix.Add("TN", scoreNX);

            matrix.Add("CA", -mismatchPenaltyBaseSubstitution);
            matrix.Add("CT", -mismatchPenaltyBaseSubstitution);
            matrix.Add("CC", scoreHit);
            matrix.Add("CG", -mismatchPenaltyBaseSubstitution);
            matrix.Add("CN", scoreNX);

            matrix.Add("GA", -mismatchPenaltyBaseSubstitution);
            matrix.Add("GT", -mismatchPenaltyBaseSubstitution);
            matrix.Add("GC", -mismatchPenaltyBaseSubstitution);
            matrix.Add("GG", scoreHit);
            matrix.Add("GN", scoreNX);

            matrix.Add("NA", scoreNX);
            matrix.Add("NT", scoreNX);
            matrix.Add("NC", scoreNX);
            matrix.Add("NG", scoreNX);
            matrix.Add("NN", scoreNX);
            return new SubstitutionMatrix(matrix, false, gapExistPenalty, gapExtendPenalty, scoreHit, -mismatchPenaltyBaseSubstitution);

        }


        /// <summary>
        /// Get a substitution matrix for nucleotide sequences.
        /// Pam25: 25% difference between the two sequences is assumed. Transitions are three times more frequent than transversions.
        /// Default values for gap exist and extend penalties are assumed: GapExistPenalty=11; GapExtendPenalty=1
        /// </summary>
        /// <returns>the substitution matrix which can be used for nucleotide sequences</returns>
        public static SubstitutionMatrix GetNucleotideSequenceSubstitutionMatrix_PAM25WithDifferentlyScoredTransitionTransversion()
        {
            return GetNucleotideSequenceSubstitutionMatrix_PAM25WithDifferentlyScoredTransitionTransversion(11.0F, 1.0F);
        }

        public static SubstitutionMatrix GetNucleotideSequenceSubstitutionMatrix_PAM25WithDifferentlyScoredTransitionTransversion(float gapExistencePenalty, float gapExtendPenalty)
        {
            float hit = 1.66F;
            float ts = -1.06F;
            float tv = -2.46F;


            Dictionary<string, float> matrix = new Dictionary<string, float>();
            matrix.Add("AA", hit);
            matrix.Add("AT", tv);
            matrix.Add("AC", tv);
            matrix.Add("AG", ts);
            matrix.Add("AN", scoreNX);

            matrix.Add("TA", tv);
            matrix.Add("TT", hit);
            matrix.Add("TC", ts);
            matrix.Add("TG", tv);
            matrix.Add("TN", scoreNX);

            matrix.Add("CA", tv);
            matrix.Add("CT", ts);
            matrix.Add("CC", hit);
            matrix.Add("CG", tv);
            matrix.Add("CN", scoreNX);

            matrix.Add("GA", ts);
            matrix.Add("GT", tv);
            matrix.Add("GC", tv);
            matrix.Add("GG", hit);
            matrix.Add("GN", scoreNX);

            matrix.Add("NA", scoreNX);
            matrix.Add("NT", scoreNX);
            matrix.Add("NC", scoreNX);
            matrix.Add("NG", scoreNX);
            matrix.Add("NN", scoreNX);
            return new SubstitutionMatrix(matrix, false, gapExistencePenalty, gapExtendPenalty, hit, tv);
        }

        public static SubstitutionMatrix GetNucleotideSequenceSubstitutionMatrix_PAM10EquallyScoredTransitionsTransversions(float gapExistencePenalty, float gapExtendPenalty)
        {
            float hit = 1.86F;
            float ts = -3.0F;
            float tv = -3.0F;


            Dictionary<string, float> matrix = new Dictionary<string, float>();
            matrix.Add("AA", hit);
            matrix.Add("AT", tv);
            matrix.Add("AC", tv);
            matrix.Add("AG", ts);
            matrix.Add("AN", scoreNX);

            matrix.Add("TA", tv);
            matrix.Add("TT", hit);
            matrix.Add("TC", ts);
            matrix.Add("TG", tv);
            matrix.Add("TN", scoreNX);

            matrix.Add("CA", tv);
            matrix.Add("CT", ts);
            matrix.Add("CC", hit);
            matrix.Add("CG", tv);
            matrix.Add("CN", scoreNX);

            matrix.Add("GA", ts);
            matrix.Add("GT", tv);
            matrix.Add("GC", tv);
            matrix.Add("GG", hit);
            matrix.Add("GN", scoreNX);

            matrix.Add("NA", scoreNX);
            matrix.Add("NT", scoreNX);
            matrix.Add("NC", scoreNX);
            matrix.Add("NG", scoreNX);
            matrix.Add("NN", scoreNX);
            return new SubstitutionMatrix(matrix, false, gapExistencePenalty, gapExtendPenalty, hit, tv);
        }

        /// <summary>
        /// Get a PAM10 substitution matrix for nucleotide sequences. 
        /// Transitions and transversions are equally scored.
        /// The gap_exist_penalty is per default 11, and the gap_extend_penalty 2.
        /// </summary>
        /// <returns></returns>
        public static SubstitutionMatrix GetNucleotideSequenceSubstitutionMatrix_PAM10EquallyScoredTransitionsTransversions()
        {
            return GetNucleotideSequenceSubstitutionMatrix_PAM10EquallyScoredTransitionsTransversions(11.0F, 2.0F);
        }

    }

        /// <summary>
        /// Class for the construction of a blastn two sequence alignment
        /// sequence 1: database
        /// sequence 2: query
        /// insertions are shifts(-) in the database-sequence, because this indicates that a insertion occured in the query with respect to the database sequence
        /// deletions are shifts(-) in the query-sequence, since this indicates that a base is missing with respect to the database sequence
        /// </summary>
        public class PairwiseAlignmentBuilder
        {

            //Downstream-> 3' direction all polymerases are moving in this direction
            //upstream->5'
            private PairwiseAlignment coreAlignment = new PairwiseAlignment(SequenceFactory.GetDefaultSequence(), SequenceFactory.GetDefaultSequence());
            private Stack<byte> databaseSequence_5_prime = new Stack<byte>();
            private Stack<byte> querySequence_5_prime = new Stack<byte>();
            private Queue<byte> databaseSequence_3_prime = new Queue<byte>();
            private Queue<byte> querySequence_3_prime = new Queue<byte>();


            public PairwiseAlignmentBuilder()
            {

            }
            public PairwiseAlignmentBuilder(PairwiseAlignment startAlignment)
            {
                this.coreAlignment = startAlignment;
            }

            /// <summary>
            /// Append a database and a query character 3' (downstream) to the sequence initiation site (eg. BLAST-seed)
            /// (uses a generic Queue to store the characters)
            /// </summary>
            public void Append_3_prime(char databaseChar, char queryChar)
            {
                databaseSequence_3_prime.Enqueue((byte)databaseChar);
                querySequence_3_prime.Enqueue((byte)queryChar);
            }


            public void Remove_3_prime(int count)
            {
                //Unfortunately the easist way for a 3 prime removal is to just get a subalignment

                //First create the new core
                PairwiseAlignment newCore = this.GetAlignment();
                ISequenceContainer data = newCore.DatabaseSequence;
                ISequenceContainer query = newCore.QuerySequence;
                coreAlignment = new PairwiseAlignment(data.SubSequence(0, data.Length - count), query.SubSequence(0, query.Length - count));
            }


            public void Remove_5_prime(int count)
            {
                //First create the new core
                PairwiseAlignment newCore = this.GetAlignment();
                ISequenceContainer data = newCore.DatabaseSequence;
                ISequenceContainer query = newCore.QuerySequence;
                coreAlignment = new PairwiseAlignment(data.SubSequence(count, data.Length - count), query.SubSequence(count, query.Length - count));


            }

            /// <summary>
            /// Appends an alignment at 3'-direction (downstream to the given Alignment);
            /// </summary>
            /// <param name="toAppend"></param>
            public void Append_3_prime(PairwiseAlignment toAppend)
            {
                if (toAppend != null)
                {

                    for (int i = 0; i < toAppend.Length; i++)
                    {
                        Append_3_prime(toAppend.DatabaseSequence[i], toAppend.QuerySequence[i]);
                    }
                }
            }


            /// <summary>
            /// Appends a database and a query character 5' (upstream) to the sequence inititation site (eg. BLAST-seed)
            /// (uses a generic Stack to store the characters)
            /// </summary>
            public void Append_5_prime(char databaseChar, char queryChar)
            {
                databaseSequence_5_prime.Push((byte)databaseChar);
                querySequence_5_prime.Push((byte)queryChar);
            }

            public void Append_5_prime(PairwiseAlignment toAppend)
            {
                if (toAppend != null)
                {
                    for (int i = toAppend.Length - 1; i >= 0; i--)
                    {
                        Append_5_prime(toAppend.DatabaseSequence[i], toAppend.QuerySequence[i]);
                    }
                }
            }

            /// <summary>
            /// Retrieves a Alignment_TwoSequences from the previously appendet characters.
            /// First the 5'-characters are assembled than the 3'-characters
            /// (pop the Stack and dequeue the Queue)
            /// </summary>
            /// <returns></returns>
            public PairwiseAlignment GetAlignment()
            {
                //First test if the length of the different sequences are valid
                if (databaseSequence_5_prime.Count != querySequence_5_prime.Count ||
                    databaseSequence_3_prime.Count != querySequence_3_prime.Count
                    ) throw new InvalidOperationException("The database-sequence and the query-sequence have to have identical length, otherwise the two-sequence-alignment is invalid");

                //If the 5' and the 3' buffer are empty return the previous alignment unmodified
                if (databaseSequence_5_prime.Count == 0 && databaseSequence_3_prime.Count == 0) return coreAlignment;


                ISequenceContainer databaseSequence = SequenceFactory.GetDefaultSequence();
                ISequenceContainer querySequence = SequenceFactory.GetDefaultSequence();

                int count_5p = databaseSequence_5_prime.Count;
                int count_3p = databaseSequence_3_prime.Count;

                //Add the 5 prime alignment
                for (int i = 0; i < count_5p; i++)
                {
                    databaseSequence.Append(databaseSequence_5_prime.Pop());
                    querySequence.Append(querySequence_5_prime.Pop());
                }


                //Add the core alignment
                databaseSequence.Append(coreAlignment.DatabaseSequence);
                querySequence.Append(coreAlignment.QuerySequence);


                //Add the 3 prime alignment
                for (int i = 0; i < count_3p; i++)
                {
                    databaseSequence.Append(databaseSequence_3_prime.Dequeue());
                    querySequence.Append(querySequence_3_prime.Dequeue());
                }

                coreAlignment = new PairwiseAlignment(databaseSequence, querySequence);
                databaseSequence_3_prime = new Queue<byte>();
                databaseSequence_5_prime = new Stack<byte>();
                querySequence_3_prime = new Queue<byte>();
                querySequence_5_prime = new Stack<byte>();
                return coreAlignment;

            }

        }



        /// <summary>
        /// this interface has to be implemented by each dynamic programming algorithm such as Smith-Waterman or Needleman-Wunsch
        /// </summary>
        public interface IDynamicProgramming
        {

            PairwiseAlignment GetAlignment();
            float Score { get;}
            SubstitutionMatrix SubstitutionMatrix { get;}

            int Start_Database { get;}
            int Start_Query { get;}
            int End_Database { get;}
            int End_Query { get;}
        }

        /// <summary>
        /// this interface has to be implemented by each pair of similar sequences,
        /// </summary>
        public interface IPairwiseAlignmentContainer : IFeature
        {
            PairwiseAlignment Alignment { get;}
            SubstitutionMatrix SubstitutionMatrix { get;set;}
            Type AlgorithmUsedForAlignment { get; set;}

            //Main info
            string DatabaseParent { get;set;}
            string QueryParent { get;set;}
            bool PlusPlusStrand { get;set;}
            int LengthAlignedDatabase { get;}
            int LengthAlignedQuery { get;}
            int? LengthDatabaseParent { get;set;}
            int? LengthQueryParent { get;set;}
            int StartDatabase { get;}
            int StartQuery { get;}
            int EndDatabase { get;}
            int EndQuery { get;}
            object Tag { get;set;}


            //Similarity gaps etc
            int Gaps { get;}
            float? Score { get;set;}
            float SimilarityWithoutGaps { get;}
            float SimilarityWithGaps { get;}
            int LengthAlignmentWithGaps { get;}
            int LengthAlignmentWithoutGaps { get;}

            //For composite
            IPairwiseAlignmentContainer this[int index] { get;}
            int Count_SubAlignments { get;}
            int Count_LongGaps(int minGapLength);
            int Count_LongGaps_Database(int minGapLength);
            int Count_LongGaps_Query(int minGapLength);
            List<IPairwiseAlignmentContainer> SubAlignments { get;}
            List<GapInPairwiseNucleotideSequenceAlignment> LongGaps(int minGapLength);
            List<GapInPairwiseNucleotideSequenceAlignment> LongGaps_Database(int minGapLength);
            List<GapInPairwiseNucleotideSequenceAlignment> LongGaps_Query(int minGapLength);

            Bio.Blast.ISignificator Significator { get;set;}
            bool IsSignificant();

            //Dubious
            int Count_TotalQueryHits { get;set;}
        }

        public interface IPairwiseNucleotideSequenceAlignmentContainer : IPairwiseAlignmentContainer
        {
            int Transitions { get;}
            int Transversions { get;}

            /// <summary>
            /// Get or set a SNPExtractor;
            /// May also be used for SNPAllele extraktion since this strategies implement the same interface
            /// </summary>
            //Bio.SNP.ISNPExtractor SNPExtractor { get;set;}
            //List<Bio.SNP.SNP> GetSNPs();
            PairwiseNucleotideSequenceAlignment GetAlignmentCoveringDatabasePosition(int databasePosition);
            PairwiseNucleotideSequenceAlignment GetAlignmentCoveringQueryPosition(int queryPosition);
        }


        /// <summary>
        /// Implementation of the Smith Waterman algorithm, modified according to Gotoh
        /// This implementation (Gotoh) is much faster than the original Smith Waterman (M*M*N) and requires only a fraction of the memory (M*N)
        /// Returns only the highest scoring pair of sequences, the second and third best can not be obtained
        /// uses the equation: Wk= PExi - PExt * (x-1) for calculating the gap penalties
        /// </summary>
        public class SmithWatermanGotoh : IDynamicProgramming
        {
            protected enum LastMove : byte { Diagonal = 1, Insertion = 2, Deletion = 3, None = 4 };
            protected struct MatrixPosition
            {
                public float score;
                public LastMove lastMove;
                public MatrixPosition(float score, LastMove lastMove)
                {
                    this.score = score;
                    this.lastMove = lastMove;
                }
            }


            //
            //CORE Variables algorithm
            //
            //The two dimensional matrix holding the best score for each position
            protected MatrixPosition[,] matrix = null;
            //Gotoh implementation additionally requires two one dimensional arrays, holding the higest deletion or insertion score
            protected float[] Dk_1; //highest score for database sequence -> insertion
            protected float[] Qk_1; //highest score for query sequence -> deletion

            //
            //Core Variables
            //
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


            //Gap penalties and default values
            protected float penaltyGapExist;
            protected float penaltyGapExtend;
            protected SmithWatermanGotoh()
            {
            }

            public SmithWatermanGotoh(ISequenceContainer databaseSequence, ISequenceContainer querySequence, SubstitutionMatrix substitutionMatrix)
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
            protected virtual void InitializeMatrix()
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
                //Call the traceback method
                Traceback();
            }
            protected float Max(float first, float second)
            {
                if (first > second) return first;
                else return second;
            }

            protected virtual void Traceback()
            {

                PairwiseAlignmentBuilder builder = new PairwiseAlignmentBuilder();
                //Actual position in the matrix; D..database; Q..query
                int posD = this.endAlignmentDatabase;
                int posQ = this.endAlignmentQuery;

                while (matrix[posD, posQ].lastMove != LastMove.None)
                {

                    //Move to the next character
                    switch (matrix[posD, posQ].lastMove)
                    {
                        case LastMove.Diagonal: //while
                            builder.Append_5_prime(databaseSequence[posD - 1], querySequence[posQ - 1]);
                            posD--; posQ--;

                            break;
                        case LastMove.Deletion:
                            int startPosD = posD;
                            do
                            {
                                builder.Append_5_prime(databaseSequence[posD - 1], '-');
                                posD--;
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


            }

            /// <summary>
            /// Get the optimal local alignment of the two specified sequences
            /// </summary>
            /// <returns></returns>
            public PairwiseAlignment GetAlignment()
            {

                return alignment;

            }


            protected static MatrixPosition GetMaximumPosition(float scoreDiagonal, float scoreInsertion, float scoreDeletion, float scoreNone)
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

                            position = new MatrixPosition(scoreDiagonal, LastMove.Diagonal);
                        }
                        else
                        {
                            //exclude scoreNone & scoreInsertion & scoreDiagonal => DELETION
                            position = new MatrixPosition(scoreDeletion, LastMove.Deletion);

                        }
                    }
                    else
                    {
                        //exclude scoreNone & scoreDiagonal


                        if (scoreInsertion > scoreDeletion)
                        {
                            //exclude scoreNone & scoreDiagonal & scoreDeletion => INSERTION
                            position = new MatrixPosition(scoreInsertion, LastMove.Insertion);
                        }
                        else
                        {
                            //exclude scoreNone &scoreDiagonal & scoreInsertion => DELETION
                            position = new MatrixPosition(scoreDeletion, LastMove.Deletion);

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
                            position = new MatrixPosition(scoreInsertion, LastMove.Insertion);
                        }
                        else
                        {
                            //exclude scoreDiagonal & scoreNone & scoreInsertion => DELETION
                            position = new MatrixPosition(scoreDeletion, LastMove.Deletion);

                        }
                    }
                    else
                    {
                        //exclude scoreDiagonal & scoreInsertion
                        if (scoreDeletion > scoreNone)
                        {
                            //exclude scoreDiagonal & scoreInsertion & scoreNone => DELETION
                            position = new MatrixPosition(scoreDeletion, LastMove.Deletion);

                        }
                        else
                        {
                            //exclude scoreDiagonal & scoreInsertion & scoreDeletion =>NONE
                            position = new MatrixPosition(scoreNone, LastMove.None);
                        }
                    }


                }


                return position; //That was annoying
            }


            protected static MatrixPosition GetMaximumPosition(float scoreDiagonal, float scoreInsertion, float scoreDeletion)
            {
                MatrixPosition position;


                //exclude scoreNone

                if (scoreDiagonal >= scoreInsertion)
                {
                    //exclude scoreNone & scoreInsertion

                    if (scoreDiagonal >= scoreDeletion)
                    {
                        //exclude scoreNone & scoreInsertion & scoreDeletion => DIAGONAL

                        position = new MatrixPosition(scoreDiagonal, LastMove.Diagonal);
                    }
                    else
                    {
                        //exclude scoreNone & scoreInsertion & scoreDiagonal => DELETION
                        position = new MatrixPosition(scoreDeletion, LastMove.Deletion);

                    }
                }
                else
                {
                    //exclude scoreNone & scoreDiagonal


                    if (scoreInsertion > scoreDeletion)
                    {
                        //exclude scoreNone & scoreDiagonal & scoreDeletion => INSERTION
                        position = new MatrixPosition(scoreInsertion, LastMove.Insertion);
                    }
                    else
                    {
                        //exclude scoreNone &scoreDiagonal & scoreInsertion => DELETION
                        position = new MatrixPosition(scoreDeletion, LastMove.Deletion);

                    }
                }



                return position; //That was annoying
            }

            /// <summary>
            /// Returns the score of the best alignment found for the two specified sequences
            /// </summary>
            public virtual float Score
            {
                get
                {
                    return highScore;
                }
            }


            public virtual int Start_Database
            {
                get
                {
                    return this.startAlignmentDatabase;
                }
            }

            /// <summary>
            /// Start postion of the optimal alignment with respect to the query sequence
            /// </summary>
            public virtual int Start_Query
            {
                get
                {
                    return this.startAlignmentQuery;
                }
            }

            /// <summary>
            /// End position of the optimal alignment with respect to the database sequence
            /// </summary>
            public virtual int End_Database
            {
                get
                {
                    return this.endAlignmentDatabase;
                }
            }

            /// <summary>
            /// End position of the optimal alignment with respect to the query sequence
            /// </summary>
            public virtual int End_Query
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
        /// Implementation of the Smith Waterman algorithm, modified according to Gotoh.
        /// Specialy adapted for the 454 sequencing, uses boundary cross penalties for gaps extending beyond poly-N tracts
        /// The direction of the 454 read is important, the query sequence has to have the 5'->3' direction 
        /// P..Plus/Plus
        /// </summary>
        public class SmithWatermanGotoh_454P : SmithWatermanGotoh
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



            protected SmithWatermanGotoh_454P()
            {
            }

            public SmithWatermanGotoh_454P(ISequenceContainer databaseSequence, ISequenceContainer querySequence, SubstitutionMatrix substitutionMatrix)
            {
                if (databaseSequence.Length < 1 || querySequence.Length < 1) throw new InvalidSequenceException("Error in Smith-Waterman-Gotoh dynamic programming. Length of database sequence or query sequence has to be at least 1 bp");


                this.penaltyGapExist = substitutionMatrix.GapExistPenalty;
                this.penaltyGapExtend = substitutionMatrix.GapExtendPenalty;
                this.databaseSequence = databaseSequence;
                this.querySequence = querySequence;
                this.substitutionMatrix = substitutionMatrix;
                this.boundaryCrossPenalty = 2 * penaltyGapExtend;

                InitializeMatrix();
            }

            public SmithWatermanGotoh_454P(ISequenceContainer databaseSequence, ISequenceContainer querySequence, SubstitutionMatrix substitutionMatrix, float homopolymereTransgressionPenalty)
            {
                if (databaseSequence.Length < 1 || querySequence.Length < 1) throw new InvalidSequenceException("Error in Smith-Waterman-Gotoh dynamic programming. Length of database sequence or query sequence has to be at least 1 bp");


                this.penaltyGapExist = substitutionMatrix.GapExistPenalty;
                this.penaltyGapExtend = substitutionMatrix.GapExtendPenalty;
                this.databaseSequence = databaseSequence;
                this.querySequence = querySequence;
                this.substitutionMatrix = substitutionMatrix;
                this.boundaryCrossPenalty = homopolymereTransgressionPenalty;

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

                while (matrix[posD, posQ].lastMove != LastMove.None)
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

                if (endAlignmentDatabase > 0 && endAlignmentQuery > 0)
                {
                    startAlignmentDatabase = posD + 1;
                    startAlignmentQuery = posQ + 1;
                    alignment = builder.GetAlignment();
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
        public class SmithWatermanGotoh_454M : SmithWatermanGotoh
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



            protected SmithWatermanGotoh_454M()
            {
            }

            public SmithWatermanGotoh_454M(ISequenceContainer databaseSequence, ISequenceContainer querySequence, SubstitutionMatrix substitutionMatrix)
            {
                if (databaseSequence.Length < 1 || querySequence.Length < 1) throw new InvalidSequenceException("Error in Smith-Waterman-Gotoh dynamic programming. Length of database sequence or query sequence has to be at least 1 bp");


                this.penaltyGapExist = substitutionMatrix.GapExistPenalty;
                this.penaltyGapExtend = substitutionMatrix.GapExtendPenalty;
                this.databaseSequence = databaseSequence;
                this.querySequence = querySequence;
                this.substitutionMatrix = substitutionMatrix;
                this.boundaryCrossPenalty = 2 * penaltyGapExtend;



                InitializeMatrix();
            }

            public SmithWatermanGotoh_454M(ISequenceContainer databaseSequence, ISequenceContainer querySequence, SubstitutionMatrix substitutionMatrix, float homopolymereTransgressionPenalty)
            {
                if (databaseSequence.Length < 1 || querySequence.Length < 1) throw new InvalidSequenceException("Error in Smith-Waterman-Gotoh dynamic programming. Length of database sequence or query sequence has to be at least 1 bp");


                this.penaltyGapExist = substitutionMatrix.GapExistPenalty;
                this.penaltyGapExtend = substitutionMatrix.GapExtendPenalty;
                this.databaseSequence = databaseSequence;
                this.querySequence = querySequence;
                this.substitutionMatrix = substitutionMatrix;
                this.boundaryCrossPenalty = homopolymereTransgressionPenalty;



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
                    matrix[i, querySequence.Length] = new MatrixPosition(0, LastMove.None);
                    Dk_1[i] = 0 - penaltyGapExist;
                    Dk_1_var[i] = 0;
                }


                //Second set all values in the vertical (query) to zero
                for (int i = 0; i < matrix.GetLength(1); i++)
                {
                    matrix[databaseSequence.Length, i] = new MatrixPosition(0, LastMove.None);
                    Qk_1[i] = 0 - penaltyGapExist;
                    Qk_1_var[i] = 0;
                }

                ProcessGapExistPenalty(); //Create position specific gap existence penalty matrix

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








                        matrix[i, k] = GetMaximumPosition(scoreDiagonal, downScoreInsertion, rightScoreDeletion, 0);

#if DEBUG
                        MatrixPosition mp = GetMaximumPosition(scoreDiagonal, downScoreInsertion, rightScoreDeletion, 0);
#endif

                        //Updating the highest scoring matrix entry
                        if (matrix[i, k].score > highScore)
                        {
                            //new highscore
                            highScore = matrix[i, k].score;
                            startAlignmentDatabase = i;
                            startAlignmentQuery = k;
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
                        float step = (penaltyGapExist - lowestGapExist) / (float)(k - 1);

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
                        float step = (penaltyGapExist - lowestGapExist) / (float)(k - 1);

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
                int posD = this.startAlignmentDatabase;
                int posQ = this.startAlignmentQuery;
                float gp = penaltyGapExist;

                while (matrix[posD, posQ].lastMove != LastMove.None)
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
                if (highScore > 0)
                {
                    endAlignmentDatabase = posD;
                    endAlignmentQuery = posQ;
                    alignment = builder.GetAlignment();
                    startAlignmentQuery++;
                    startAlignmentDatabase++;
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
        /// Implementation of the Needleman Wunsch algorithm, modified according to Gotoh;
        /// uses the equation: Wk= PExi - PExt * (x-1) for calculating the gap penalties
        /// </summary>
        public class NeedlemanWunschGotoh : SmithWatermanGotoh
        {


            public NeedlemanWunschGotoh(ISequenceContainer databaseSequence, ISequenceContainer querySequence, SubstitutionMatrix substitutionMatrix)
            {
                if (databaseSequence.Length < 1 || querySequence.Length < 1) throw new InvalidSequenceException("Error in Needleman-Wunsch-Gotoh dynamic programming. Length of database sequence or query sequence has to be at least 1 bp");
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
                this.Qk_1 = new float[querySequence.Length + 1];

                //First the zeroPosition
                matrix[0, 0] = new MatrixPosition(0, LastMove.Diagonal);
                Dk_1[0] = 0.0F;
                Qk_1[0] = 0.0F;


                //Second set all values in the horizontal (database) to the initate value minus the penalties
                for (int i = 1; i < matrix.GetLength(0); i++)
                {
                    float valueAttempt = 0.0F - penaltyGapExist - (penaltyGapExtend * (i - 1));

                    matrix[i, 0] = new MatrixPosition(valueAttempt, LastMove.Deletion);
                    //leafe as zero
                    Dk_1[i] = valueAttempt - penaltyGapExist;
                }

                //Third set all values in the vertical (query) to zero
                for (int i = 1; i < matrix.GetLength(1); i++)
                {
                    float valueAttempt = 0.0F - penaltyGapExist - (penaltyGapExtend * (i - 1));


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



                    }
                }
                Traceback();
            }

            protected override void Traceback()
            {
                //If the matrix has not yet been constructed, initialize it
                if (matrix == null) InitializeMatrix();

                PairwiseAlignmentBuilder builder = new PairwiseAlignmentBuilder();
                //Actual position in the matrix; D..database; Q..query
                int posD = databaseSequence.Length;
                int posQ = querySequence.Length;

                while (posD != 0 || posQ != 0)
                {

                    //Move to the next character
                    switch (matrix[posD, posQ].lastMove)
                    {
                        case LastMove.Diagonal: //while
                            builder.Append_5_prime(databaseSequence[posD - 1], querySequence[posQ - 1]);
                            posD--; posQ--;

                            break;
                        case LastMove.Deletion:
                            int startPosD = posD;
                            do
                            {
                                builder.Append_5_prime(databaseSequence[posD - 1], '-');
                                posD--;
#if DEBUG
                                double a = Math.Round(matrix[posD, posQ].score, 3);
                                double b = Math.Round((matrix[startPosD, posQ].score + penaltyGapExist + penaltyGapExtend * (startPosD - posD - 1)), 3);
#endif

                            }
                            while (Math.Round(matrix[posD, posQ].score, 3) != Math.Round((matrix[startPosD, posQ].score + penaltyGapExist + penaltyGapExtend * (startPosD - posD - 1)), 3));

                            break;
                        case LastMove.Insertion:
                            int startPosQ = posQ;
                            do
                            {
                                builder.Append_5_prime('-', querySequence[posQ - 1]);
                                posQ--;

                            }
                            while (Math.Round(matrix[posD, posQ].score, 2) != Math.Round((matrix[posD, startPosQ].score + penaltyGapExist + penaltyGapExtend * (startPosQ - posQ - 1)), 2));


                            break;
                        default: throw new InvalidOperationException("this should not happen");
                    }


                }



                startAlignmentDatabase = posD + 1;
                startAlignmentQuery = posQ + 1;

                this.endAlignmentDatabase = databaseSequence.Length;
                this.endAlignmentQuery = querySequence.Length;
                this.highScore = matrix[endAlignmentDatabase, endAlignmentQuery].score;

                alignment = builder.GetAlignment();

            }






        }

        public class PairwiseNucleotideSequenceAlignment : IPositionable, IPairwiseNucleotideSequenceAlignmentContainer
        {


            private PairwiseAlignment alignment = null;
            private string databaseParentName = null;
            private string queryParentName = null;
            private int startPosDatabase = 0;
            private int endPosDatabase = 0;
            private int startPosQuery = 0;
            private int endPosQuery = 0;
            private bool? strandPlusPlus = null;
            private int? lengthDatabaseParent = null;
            private object tag;
            private int? lengthQueryParent = null;


            private Bio.Blast.ISignificator significator = null;

            /// <summary>
            /// Settings
            /// </summary>
            private SubstitutionMatrix substitutionMatrix = null;
            private Type algorithmUsedForAlignment = null;
            private float? score;




            //WORKING VARIABLES
            //Transitions A-G T-C
            //Transversions A-C A-T C-G 
            private int? hits = null;
            private int? baseSubstitutions = null;
            private float? similarity_withoutGaps = null;
            private float? similarity_withGaps = null;
            private int? transitions = null;
            private int? transversions = null;
            private int? indels = null;
            private int? insertions = null;
            private int? deletions = null;
            private int? lengthDatabase = null;
            private int? lengthQuery = null;
            private int? lengthAlignmentWithoutGaps = null;
            private int? gaps = null;

            //Dubious
            private int countTotalQueryHits = 1;


            public PairwiseNucleotideSequenceAlignment(PairwiseAlignment alignment, string databaseParentName, string queryParentName, int startPostionDatabase, int startPositionQuery, int endPositionDatabase, int endPositionQuery)
            {
                this.alignment = alignment;
                this.databaseParentName = databaseParentName;
                this.queryParentName = queryParentName;
                this.startPosDatabase = startPostionDatabase;
                this.startPosQuery = startPositionQuery;
                this.endPosDatabase = endPositionDatabase;
                this.endPosQuery = endPositionQuery;
            }

            public PairwiseNucleotideSequenceAlignment(PairwiseAlignment alignment, string databaseParentName, string queryParentName, int startPostionDatabase, int startPositionQuery, int endPositionDatabase, int endPositionQuery, bool plusPlusStrand)
            {
                this.alignment = alignment;
                this.databaseParentName = databaseParentName;
                this.queryParentName = queryParentName;
                this.startPosDatabase = startPostionDatabase;
                this.startPosQuery = startPositionQuery;
                this.endPosDatabase = endPositionDatabase;
                this.endPosQuery = endPositionQuery;
                this.strandPlusPlus = plusPlusStrand;
            }




            /// <summary>
            /// The score achieved by the alignment
            /// </summary>
            public float? Score
            {
                get
                {
                    return score;
                }
                set
                {
                    this.score = value;
                }
            }


            /// <summary>
            /// The total number of gaps in the sequence, irrespective whether the gap is in the database or in the query sequence.
            /// Actually is the number of occurences of the char '-'
            /// </summary>
            public int Gaps
            {
                get
                {
                    if (this.gaps == null)
                    {
                        this.gaps = 0;
                        ISequenceContainer data = alignment.DatabaseSequence;
                        ISequenceContainer query = alignment.QuerySequence;

                        for (int i = 0; i < data.Length; i++)
                        {
                            if (data[i] == '-' || query[i] == '-') gaps++;
                        }
                    }
                    return this.gaps.Value;
                }

            }





            /// <summary>
            /// Is the start position of the alignment in the database sequence, value is identical with StartDatabase property
            /// </summary>
            public int? Start
            {
                get
                {
                    return this.StartDatabase;
                }
                set
                {
                    throw new InvalidOperationException("Setting of Start position is not valid for alignment; Start positions are obtained through dynamic programming");
                }
            }

            /// <summary>
            /// Same as EndDatabase property; 
            /// setting is not permitted 
            /// </summary>
            /// <exception cref="InvalidOperationException">when attempting to set this property to a certain value an exception will be thrown</exception>
            public int? End
            {
                get
                {
                    return this.EndDatabase;
                }
                set
                {
                    throw new InvalidOperationException("Setting of End position is not valid for alignments. End positions are obtained through dynamic programming");
                }
            }
            /// <summary>
            /// Is this sequence the root? For alignments always false; Set property throws error
            /// </summary>
            /// <exception cref="InvalidOperationException">when attempting to set this property to a certain value an exception will be thrown</exception>
            public bool? IsRoot
            {
                get
                {
                    return false;
                }
                set
                {
                    throw new InvalidOperationException("Setting of IsRoot property is not valid for alignments; An alignment can never be the root");
                }
            }

            /// <summary>
            /// Same as DatabaseParent property
            /// </summary>
            public string ParentName
            {
                get
                {

                    return this.DatabaseParent;
                }
                set
                {
                    this.DatabaseParent = value;
                }
            }

            /// <summary>
            /// Same as LengthDatabase property
            /// </summary>
            public long Length
            {
                get
                {
                    return this.LengthAlignedDatabase;
                }
            }


            /// <summary>
            /// Retrieve a Subalignment. From the start position to the end of the alignment, relative to the query sequence 
            /// </summary>
            /// <param name="startPostion">start position of the subalignment with respect to the query sequence</param>
            /// <returns></returns>
            public PairwiseNucleotideSequenceAlignment SubAlignment_RelativeToQuery(int startPostion)
            {
                PairwiseAlignmentBuilder builder = new PairwiseAlignmentBuilder();
                PairwiseAlignmentBuilder buildCleav = new PairwiseAlignmentBuilder();
                ISequenceContainer query = alignment.QuerySequence;
                ISequenceContainer database = alignment.DatabaseSequence;

                int realRelativeStartQuery = startPostion - startPosQuery;
                int realRelativeStartDatabase = 0;
                int relativeEndPosDatabase = 0;
                int relativeEndPosQuery = 0;

                int countQuery = -1;
                int countDatabase = -1;



                //Calculate real relative StartPosition
                bool append = false;
                for (int i = 0; i < alignment.Length; i++)
                {
                    //Increase the relative counts of the query and the database sequence if the character is not a gap
                    if (query[i] != '-') countQuery++;
                    if (database[i] != '-') countDatabase++;

                    //has the start of the subalignment been reached
                    if (!append && countQuery >= realRelativeStartQuery && database[i] != '-' && query[i] != '-')
                    {
                        append = true;
                        realRelativeStartQuery = countQuery;
                        realRelativeStartDatabase = countDatabase;

                    }

                    /*/has the end of the subalignment been reached
                    if (realRelativeStartQuery + Length == countQuery)
                    {
                        append = false;
                        relativeEndPosQuery = countQuery;
                        relativeEndPosDatabase = countDatabase;
                        break;

                    }*/

                    //Create the alignment which will be cleaved
                    if (!append)
                    {
                        buildCleav.Append_3_prime(database[i], query[i]);
                    }

                    if (append)
                    {
                        builder.Append_3_prime(database[i], query[i]);

                    }
                }


                PairwiseAlignment align = builder.GetAlignment();

                PairwiseNucleotideSequenceAlignment al = new PairwiseNucleotideSequenceAlignment(align, this.DatabaseParent, this.QueryParent
                    , this.StartDatabase + realRelativeStartDatabase
                    , this.StartQuery + realRelativeStartQuery
                    , this.EndDatabase
                    , this.EndQuery);
                al.PlusPlusStrand = this.PlusPlusStrand;
                al.SubstitutionMatrix = this.SubstitutionMatrix;
                al.AlgorithmUsedForAlignment = this.AlgorithmUsedForAlignment;
                if (al.SubstitutionMatrix != null) al.Score = this.Score - al.SubstitutionMatrix.GetScoreForAlignment(buildCleav.GetAlignment());
                else al.Score = this.Score;
                return al;


            }

            /// <summary>
            /// Retrieve a Subalignment. Start position and length are relative to the query sequence
            /// </summary>
            /// <param name="startPostion">start position of the subalignment with respect to the query sequence</param>
            /// <returns></returns>
            public PairwiseNucleotideSequenceAlignment SubAlignment_RelativeToQuery(int startPostion, int length)
            {
                PairwiseAlignmentBuilder builder = new PairwiseAlignmentBuilder();
                PairwiseAlignmentBuilder buildCleav = new PairwiseAlignmentBuilder();
                ISequenceContainer query = alignment.QuerySequence;
                ISequenceContainer database = alignment.DatabaseSequence;

                int realRelativeStartQuery = startPostion - startPosQuery;
                int realRelativeStartDatabase = 0;
                int relativeEndPosDatabase = 0;
                int relativeEndPosQuery = 0;

                int countQuery = -1;
                int countDatabase = -1;



                //Calculate real relative StartPosition
                bool append = false;
                bool end = false;
                for (int i = 0; i < alignment.Length; i++)
                {
                    //Increase the relative counts of the query and the database sequence if the character is not a gap
                    if (query[i] != '-') countQuery++;
                    if (database[i] != '-') countDatabase++;

                    //has the start of the subalignment been reached
                    if (!append && !end && countQuery >= realRelativeStartQuery && database[i] != '-' && query[i] != '-')
                    {
                        append = true;
                        realRelativeStartQuery = countQuery;
                        realRelativeStartDatabase = countDatabase;

                    }

                    //has the end of the subalignment been reached
                    if (realRelativeStartQuery + Length == countQuery)
                    {
                        end = true;
                        relativeEndPosQuery = countQuery;
                        relativeEndPosDatabase = countDatabase;
                        break;

                    }

                    //Create the alignment which will be cleaved
                    if (!append && !end)
                    {
                        buildCleav.Append_3_prime(database[i], query[i]);
                    }

                    if (append && !end)
                    {
                        builder.Append_3_prime(database[i], query[i]);

                    }
                }


                PairwiseAlignment align = builder.GetAlignment();

                PairwiseNucleotideSequenceAlignment al = new PairwiseNucleotideSequenceAlignment(align, this.DatabaseParent, this.QueryParent
                    , this.StartDatabase + realRelativeStartDatabase
                    , this.StartQuery + realRelativeStartQuery
                    , this.StartDatabase + relativeEndPosDatabase
                    , this.StartQuery + relativeEndPosQuery);
                al.PlusPlusStrand = this.PlusPlusStrand;
                al.SubstitutionMatrix = this.SubstitutionMatrix;
                al.AlgorithmUsedForAlignment = this.AlgorithmUsedForAlignment;
                if (this.SubstitutionMatrix != null) al.Score = this.Score - al.SubstitutionMatrix.GetScoreForAlignment(buildCleav.GetAlignment());
                else al.Score = this.Score;
                return al;







            }

            /// <summary>
            /// Get start position of the database
            /// </summary>
            public int StartDatabase
            {
                get
                {
                    return this.startPosDatabase;
                }
            }

            /// <summary>
            /// Get the start position of the query sequence, i.e. position with respect to the query parent
            /// </summary>
            public int StartQuery
            {
                get
                {
                    return this.startPosQuery;
                }
            }

            /// <summary>
            /// Is the given alignment on the +/+ strand (true) or on the +/- strand (false)
            /// </summary>
            public bool PlusPlusStrand
            {
                get
                {
                    return this.strandPlusPlus.Value;
                }
                set
                {
                    this.strandPlusPlus = value;
                }
            }

            /// <summary>
            /// Get the end position of the database sequence, i.e. position with respect to the database parent
            /// </summary>
            public int EndDatabase
            {
                get
                {
                    return this.endPosDatabase;
                }
            }

            /// <summary>
            /// Get the end position of the query sequence, i.e. position with respect to the query parent
            /// </summary>
            public int EndQuery
            {
                get
                {
                    return this.endPosQuery;
                }
            }

            /// <summary>
            /// Return the PairwiseAlignment
            /// </summary>
            public PairwiseAlignment Alignment
            {
                get
                {
                    return this.alignment;
                }
            }

            /// <summary>
            /// Get or set the substitution matrix used for creating the alignments
            /// </summary>
            public SubstitutionMatrix SubstitutionMatrix
            {
                get
                {
                    return this.substitutionMatrix;
                }
                set
                {
                    this.substitutionMatrix = value;
                }
            }

            /// <summary>
            /// Get or set the name of the database parent sequence
            /// </summary>
            public string DatabaseParent
            {
                get
                {
                    return this.databaseParentName;
                }
                set
                {
                    this.databaseParentName = value;
                }
            }

            /// <summary>
            /// Get or set the name of the query parent sequence
            /// </summary>
            public string QueryParent
            {
                get
                {
                    return this.queryParentName;
                }
                set
                {
                    this.queryParentName = value;
                }
            }

            /// <summary>
            /// Get the length of the database sequence. Gaps are not considered
            /// </summary>
            public int LengthAlignedDatabase
            {
                get
                {
                    if (this.lengthDatabase == null)
                    {
                        this.lengthDatabase = this.endPosDatabase - this.startPosDatabase + 1;

                    }
                    return this.lengthDatabase.Value;
                }
            }

            /// <summary>
            /// Get the length of the query sequence. Gaps are not considered
            /// </summary>
            public int LengthAlignedQuery
            {
                get
                {
                    if (this.lengthQuery == null)
                    {
                        this.lengthQuery = this.endPosQuery - this.startPosQuery + 1;

                    }
                    return this.lengthQuery.Value;
                }
            }

            /// <summary>
            /// Absolte length of the alignment. Including gaps
            /// </summary>
            public int LengthAlignmentWithGaps
            {
                get
                {
                    return alignment.Length;
                }
            }

            /// <summary>
            /// Length of the alignment not considering gaps, i.e the sum of hits and mismatches
            /// </summary>
            public int LengthAlignmentWithoutGaps
            {
                get
                {
                    if (lengthAlignmentWithoutGaps == null)
                    {
                        ISequenceContainer data = alignment.DatabaseSequence;
                        ISequenceContainer query = alignment.QuerySequence;
                        int countValid = 0;
                        for (int i = 0; i < alignment.Length; i++)
                        {
                            if (data[i] == '-' || query[i] == '-') continue;
                            countValid++;
                        }
                        lengthAlignmentWithoutGaps = countValid;
                    }
                    return lengthAlignmentWithoutGaps.Value;
                }
            }


            /// <summary>
            /// Is the given alignment significant.
            /// Strategy pattern, uses a specified ISignificator
            /// </summary>
            /// <returns></returns>
            public bool IsSignificant()
            {
                return this.significator.IsSignificant(this);
            }

            /// <summary>
            /// Get or set the type of the algorithm used for creating this alignment
            /// </summary>
            public Type AlgorithmUsedForAlignment
            {
                get
                {
                    return this.algorithmUsedForAlignment;
                }
                set
                {
                    this.algorithmUsedForAlignment = value;
                }
            }

            /// <summary>
            /// The similarity between two sequences in percent. Gaps are not considered
            /// </summary>
            public float SimilarityWithoutGaps
            {
                get
                {
                    if (similarity_withoutGaps == null)
                    {
                        ISequenceContainer data = alignment.DatabaseSequence;
                        ISequenceContainer query = alignment.QuerySequence;

                        int countIdentity = 0;
                        for (int i = 0; i < alignment.Length; i++)
                        {
                            if (data[i] == '-' || query[i] == '-') continue;


                            if (data[i] == query[i]) countIdentity++;


                        }

                        this.similarity_withoutGaps = ((float)(100 * countIdentity)) / (float)this.LengthAlignmentWithoutGaps;
                    }
                    return this.similarity_withoutGaps.Value;
                }
            }


            /// <summary>
            /// The similarity between two sequences in percent. Gaps are considered. 
            /// Value is lower or equal to SimilarityWithoutGaps
            /// </summary>
            public float SimilarityWithGaps
            {
                get
                {
                    if (similarity_withGaps == null)
                    {


                        this.similarity_withGaps = ((float)(100 * this.Hits)) / (float)this.LengthAlignmentWithGaps;
                    }
                    return this.similarity_withGaps.Value;
                }
            }


            /// <summary>
            /// Get the number of hits, i.e identical bases
            /// </summary>
            public int Hits
            {
                get
                {
                    if (this.hits == null)
                    {
                        ISequenceContainer similarity = alignment.SimilaritySequence;
                        int countHits = 0;
                        for (int i = 0; i < alignment.Length; i++)
                        {
                            if (similarity[i] == '|') countHits++;
                        }

                        this.hits = countHits;
                    }
                    return this.hits.Value;
                }
            }


            /// <summary>
            /// Get the number of transitions
            /// </summary>
            public int Transitions
            {
                get
                {
                    if (transitions == null)
                    {
                        ISequenceContainer data = alignment.DatabaseSequence;
                        ISequenceContainer query = alignment.QuerySequence;
                        int counttransitions = 0;
                        for (int i = 0; i < alignment.Length; i++)
                        {
                            switch (data[i])
                            {
                                case 'A':
                                case 'a': if (query[i] == 'G' || query[i] == 'g') counttransitions++;
                                    break;
                                case 'T':
                                case 't': if (query[i] == 'C' || query[i] == 'c') counttransitions++;
                                    break;
                                case 'C':
                                case 'c': if (query[i] == 'T' || query[i] == 't') counttransitions++;
                                    break;
                                case 'G':
                                case 'g': if (query[i] == 'A' || query[i] == 'a') counttransitions++;
                                    break;
                                default: break;


                            }

                        }
                        this.transitions = counttransitions;
                    }
                    return transitions.Value;
                }
            }
            /// <summary>
            /// Get the number of transversions
            /// </summary>
            public int Transversions
            {
                get
                {
                    if (transversions == null)
                    {
                        ISequenceContainer data = alignment.DatabaseSequence;
                        ISequenceContainer query = alignment.QuerySequence;
                        int counttransversions = 0;
                        for (int i = 0; i < alignment.Length; i++)
                        {
                            if (data[i] == '-' || query[i] == '-' || data[i] == query[i]) continue;

                            switch (data[i])
                            {
                                case 'A':
                                case 'a': if (!(query[i] == 'G' || query[i] == 'g')) counttransversions++;
                                    break;
                                case 'T':
                                case 't': if (!(query[i] == 'C' || query[i] == 'c')) counttransversions++;
                                    break;
                                case 'C':
                                case 'c': if (!(query[i] == 'T' || query[i] == 't')) counttransversions++;
                                    break;
                                case 'G':
                                case 'g': if (!(query[i] == 'A' || query[i] == 'a')) counttransversions++;
                                    break;
                                default: break;


                            }

                        }
                        this.transversions = counttransversions;
                    }
                    return transversions.Value;
                }
            }






            /// <summary>
            /// Get the Subalignment with the specified index.
            /// PairwiseNucleotideSequenceAlignments do not contain subalignments.
            /// When applied with the value 0, this instance will be returned
            /// </summary>
            /// <param name="index">index of the subalignment</param>
            /// <returns></returns>
            public IPairwiseAlignmentContainer this[int index]
            {
                get
                {
                    if (index == 0) return this;
                    else throw new ArgumentOutOfRangeException("index out of range");

                }
            }
            /// <summary>
            /// Get the number of subalignments.
            /// PairwiseNucleotideSequenceAligments do not contain subalignments, therefore the value 1 will be returned
            /// </summary>
            public int Count_SubAlignments
            {
                get { return 1; }
            }

            /// <summary>
            /// Get or set the ISignificator instance which will be used to determine whether or not this PairwiseNucleotideSequencAlignment is significant
            /// </summary>
            public Bio.Blast.ISignificator Significator
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






            /// <summary>
            /// Count the number of long gaps both in the database and in the query sequence which exceed the specified minimum length;
            /// </summary>
            /// <param name="minGapLength">minimum length of the gap</param>
            /// <returns></returns>
            public int Count_LongGaps(int minGapLength)
            {
                return this.Count_LongGaps_Database(minGapLength) + this.Count_LongGaps_Query(minGapLength);

            }

            /// <summary>
            /// Count the number of long gaps in the database sequence which exceed the specified minimum length
            /// </summary>
            /// <param name="minGapLength">minimum length of the gap, greater or equal</param>
            /// <returns></returns>
            public int Count_LongGaps_Database(int minGapLength)
            {
                ISequenceContainer data = alignment.DatabaseSequence;

                int count = 0;

                for (int i = 0; i < data.Length; i++)
                {
                    int k = 0;
                    while (data[i + k] == '-') k++;
                    if (k >= minGapLength)
                    {
                        count++;
                        i = i + k - 1;
                    }
                }
                return count;
            }


            /// <summary>
            /// 
            /// </summary>
            /// <param name="minGapLength"></param>
            /// <returns></returns>
            public List<GapInPairwiseNucleotideSequenceAlignment> LongGaps_Database(int minGapLength)
            {
                ISequenceContainer data = alignment.DatabaseSequence;
                ISequenceContainer query = alignment.QuerySequence;

                List<GapInPairwiseNucleotideSequenceAlignment> gaps = new List<GapInPairwiseNucleotideSequenceAlignment>();

                int count_database = -1;
                int count_query = -1;



                for (int i = 0; i < data.Length; i++)
                {
                    if (query[i] != '-') count_query++;
                    if (data[i] != '-') count_database++;

                    int k = 0;
                    while (data[i + k] == '-') k++;

                    //012-----345678
                    //AAA-----GGGGG gap_length=5
                    //AAATTTTTGGGGG
                    //012345678901233

                    if (k >= minGapLength)
                    {
                        gaps.Add(new GapInPairwiseNucleotideSequenceAlignment(minGapLength, count_database, count_database + 1, count_query, count_query + k + 1, this));
                        i = i + k - 1;
                    }
                }
                return gaps;
            }


            /// <summary>
            /// Count the number of long gaps in the query sequence which exceed the specified minimum length
            /// </summary>
            /// <param name="minGapLength">minimum length of the gap, greater or equal</param>
            /// <returns></returns>
            public int Count_LongGaps_Query(int minGapLength)
            {
                ISequenceContainer query = alignment.QuerySequence;
                int count = 0;

                for (int i = 0; i < query.Length; i++)
                {
                    int k = 0;
                    while (query[i + k] == '-') k++;
                    if (k >= minGapLength)
                    {
                        count++;
                        i = i + k - 1;
                    }
                }
                return count;
            }

            /// <summary>
            /// 
            /// </summary>
            /// <param name="minGapLength"></param>
            /// <returns></returns>
            public List<GapInPairwiseNucleotideSequenceAlignment> LongGaps_Query(int minGapLength)
            {
                ISequenceContainer data = alignment.DatabaseSequence;
                ISequenceContainer query = alignment.QuerySequence;

                List<GapInPairwiseNucleotideSequenceAlignment> gaps = new List<GapInPairwiseNucleotideSequenceAlignment>();

                int count_database = -1;
                int count_query = -1;



                for (int i = 0; i < data.Length; i++)
                {
                    if (query[i] != '-') count_query++;
                    if (data[i] != '-') count_database++;

                    int k = 0;
                    while (query[i + k] == '-') k++;

                    //012-----345678
                    //AAA-----GGGGG gap_length=5
                    //AAATTTTTGGGGG
                    //012345678901233

                    if (k >= minGapLength)
                    {
                        gaps.Add(new GapInPairwiseNucleotideSequenceAlignment(minGapLength, count_database, count_database + k + 1, count_query, count_query + 1, this));
                        i = i + k - 1;
                    }
                }
                return gaps;
            }


            /// <summary>
            /// Return a collection of long gaps
            /// </summary>
            /// <param name="minGapLength"></param>
            /// <returns></returns>
            public List<GapInPairwiseNucleotideSequenceAlignment> LongGaps(int minGapLength)
            {
                List<GapInPairwiseNucleotideSequenceAlignment> gaps = new List<GapInPairwiseNucleotideSequenceAlignment>();
                gaps.AddRange(LongGaps_Database(minGapLength));
                gaps.AddRange(LongGaps_Query(minGapLength));
                return gaps;
            }





            /// <summary>
            /// Get or set the length of the database parent, eg. the length of the gene or chromosome
            /// </summary>
            public int? LengthDatabaseParent
            {
                get
                {

                    return this.lengthDatabaseParent;
                }
                set
                {
                    this.lengthDatabaseParent = value;
                }
            }

            /// <summary>
            /// Get or set the length of the query parent, eg. the length of a transcript
            /// </summary>
            public int? LengthQueryParent
            {
                get
                {
                    return this.lengthQueryParent;
                }
                set
                {
                    this.lengthQueryParent = value;
                }
            }

            public List<IPairwiseAlignmentContainer> SubAlignments
            {
                get
                {

                    List<IPairwiseAlignmentContainer> temp = new List<IPairwiseAlignmentContainer>();
                    temp.Add((IPairwiseAlignmentContainer)this);
                    return temp;
                }
            }





            public PairwiseNucleotideSequenceAlignment GetAlignmentCoveringDatabasePosition(int databasePosition)
            {
                if (this.StartDatabase <= databasePosition && this.EndDatabase >= databasePosition) return this;
                else return null;
            }

            public PairwiseNucleotideSequenceAlignment GetAlignmentCoveringQueryPosition(int queryPosition)
            {
                if (this.StartQuery <= queryPosition && this.EndQuery >= queryPosition) return this;
                else return null;
            }


            public object Tag
            {
                get
                {
                    return this.tag;
                }
                set
                {
                    this.tag = value;
                }
            }


            public int Count_TotalQueryHits
            {
                get
                {
                    return this.countTotalQueryHits;
                }
                set
                {
                    this.countTotalQueryHits = value;
                }
            }





     


  


            public string FeatureName
            {
                get { return "Pairwise nucleotide sequence alignment"; }
            }


        }

        /// <summary>
        /// CompositePairwiseNucleotideSequenceAlignment are container for >1 PairwiseNucleotideSequenceAlignments.
        /// They are necessary for alignments containing large gaps, e.g. introns
        /// </summary>
        public class CompositePairwiseNucleotideSequenceAlignment : IPairwiseNucleotideSequenceAlignmentContainer, IPositionable
        {


            private List<PairwiseNucleotideSequenceAlignment> alignments;
            private PairwiseAlignment pwAlignment = null;
            private string databaseParentName;
            private string queryParentName;
            private int startPositionDatabase;
            private int endPositionDatabase;
            private int startPositionQuery;
            private int endPositionQuery;
            private bool? plusplusOrientation;
            private object tag;

            private Bio.Blast.ISignificator significator;

            private SubstitutionMatrix substitutionMatrix;
            private Type algorithmUsedForIdentification;
            private int? lengthOfDatabase = null;
            private int? lengthOfQuery = null;
            private int? gaps = null;
            private float? similaritywithoutgaps = null;
            private float? similaritywithgaps = null;
            private float? score = null;
            private int? lengthWithGaps = null;
            private int? lengthWithoutGaps = null;
            private long? length = null;
            private int? transitions = null;
            private int? transversions = null;


            //Dubios
            private int countTotalQueryHits = 1;

            /// <summary>
            /// Create a CompositePairwiseNucleotideSequenceAlignment. This alignment is a container for >1 PairwiseNucleotideSequenceAlignments.
            /// They are necessary for alignments containing large gaps, e.g. introns
            /// </summary>
            /// <param name="nucleotideSequenceAlignments">a list of Sub-alignments which will constitute the CompositePairwiseNucleotideSequenceAlignment</param>
            public CompositePairwiseNucleotideSequenceAlignment(List<PairwiseNucleotideSequenceAlignment> nucleotideSequenceAlignments)
            {
                this.alignments = nucleotideSequenceAlignments;
                this.databaseParentName = alignments[0].DatabaseParent;
                this.queryParentName = alignments[0].QueryParent;
                this.plusplusOrientation = alignments[0].PlusPlusStrand;

                ///The name of the database and the query parent has to be identical for all sub alignments, furthermore the plus/plus strand also has to be identical
                for (int i = 1; i < alignments.Count; i++)
                {
                    if (alignments[i].DatabaseParent != databaseParentName || alignments[i].QueryParent != queryParentName) throw new InvalidOperationException("The names of the database and/or sequence parents have to be identical for all subsequences");
                    if (alignments[i].PlusPlusStrand != this.plusplusOrientation) throw new InvalidOperationException("All subsequences have to be from the same orientation, either plus/plus or plus/minus");
                }

               alignments.Sort(new Bio.Seq.Sort.SortPositionables_ParentName_StartPos_Ascending<PairwiseNucleotideSequenceAlignment>());
                this.startPositionDatabase = alignments[0].StartDatabase;
                this.startPositionQuery = alignments[0].StartQuery;
                this.endPositionDatabase = alignments[alignments.Count - 1].EndDatabase;
                this.endPositionQuery = alignments[alignments.Count - 1].EndQuery;

            }



            public PairwiseAlignment Alignment
            {
                get
                {
                    if (pwAlignment == null)
                    {
                        PairwiseAlignmentBuilder builder = new PairwiseAlignmentBuilder();
                        for (int i = 0; i < alignments.Count - 1; i++)
                        {
                            builder.Append_3_prime(alignments[i].Alignment);

                            if (i < alignments.Count - 1) //Fetzengeiler algorithmus!!
                            {

                                //concatenates individual Pairwise alignments, if the alignment consists of individual alignments
                                //the separating char N will be used. Maximum number of N's is 10!
                                int count_n_database = alignments[i + 1].StartDatabase - alignments[i].EndDatabase - 1;
                                int count_n_query = alignments[i + 1].StartQuery - alignments[i].EndQuery - 1;
                                count_n_database = count_n_database > 10 ? 10 : count_n_database;
                                count_n_query = count_n_query > 10 ? 10 : count_n_query;


                                while (count_n_query > 0 && count_n_database > 0)
                                {
                                    builder.Append_3_prime('N', 'N');
                                    count_n_database--;
                                    count_n_query--;
                                }
                                while (count_n_database > 0)
                                {
                                    builder.Append_3_prime('N', '-');
                                    count_n_database--;
                                }
                                while (count_n_query > 0)
                                {
                                    builder.Append_3_prime('-', 'N');
                                    count_n_query--;
                                }



                            }
                        }
                        builder.Append_3_prime(alignments[alignments.Count - 1].Alignment);
                        pwAlignment = builder.GetAlignment();
                    }
                    return pwAlignment;


                }
            }

            public SubstitutionMatrix SubstitutionMatrix
            {
                get
                {
                    return this.substitutionMatrix;
                }
                set
                {
                    this.substitutionMatrix = value;
                }
            }

            public Type AlgorithmUsedForAlignment
            {
                get
                {
                    return this.algorithmUsedForIdentification;
                }
                set
                {
                    this.algorithmUsedForIdentification = value;
                }
            }

            /// <summary>
            /// Get or set the name of the database parent name
            /// </summary>
            public string DatabaseParent
            {
                get
                {
                    return this.databaseParentName;
                }
                set
                {
                    this.databaseParentName = value;
                    foreach (PairwiseNucleotideSequenceAlignment pw in alignments)
                    {
                        pw.DatabaseParent = value;
                    }
                }
            }

            /// <summary>
            /// Get or set the name of the query parent sequence
            /// </summary>
            public string QueryParent
            {
                get
                {
                    return this.queryParentName;
                }
                set
                {
                    this.queryParentName = value;
                    foreach (PairwiseNucleotideSequenceAlignment pw in alignments)
                    {
                        pw.QueryParent = value;
                    }

                }
            }

            /// <summary>
            /// Get or Set the orientation of the alignment either +/+ (true) or +/- (false)
            /// </summary>
            public bool PlusPlusStrand
            {
                get
                {

                    return this.plusplusOrientation.Value;
                }
                set
                {
                    this.plusplusOrientation = value;
                    foreach (PairwiseNucleotideSequenceAlignment pw in alignments)
                    {
                        pw.PlusPlusStrand = value;
                    }

                }
            }

            /// <summary>
            /// Get the length of the Database sequence, gaps are not considered
            /// </summary>
            public int LengthAlignedDatabase
            {
                get
                {
                    if (this.lengthOfDatabase == null)
                    {
                        this.lengthOfDatabase = 0;
                        for (int i = 0; i < alignments.Count; i++)
                        {
                            lengthOfDatabase += alignments[i].LengthAlignedDatabase;
                        }
                    }
                    return lengthOfDatabase.Value;
                }
            }

            /// <summary>
            /// Get the length of the query sequence. Gaps are not considered
            /// </summary>
            public int LengthAlignedQuery
            {
                get
                {
                    if (this.lengthOfQuery == null)
                    {
                        this.lengthOfQuery = 0;
                        for (int i = 0; i < alignments.Count; i++)
                        {
                            lengthOfQuery += alignments[i].LengthAlignedQuery;
                        }
                    }
                    return lengthOfQuery.Value;
                }
            }

            /// <summary>
            /// Get the start postion of the database sequence, with respect to the database parent
            /// </summary>
            public int StartDatabase
            {
                get
                {
                    return this.startPositionDatabase;
                }
            }
            /// <summary>
            /// Get the start position of the query sequence with respect to the query parent
            /// </summary>
            public int StartQuery
            {
                get
                {
                    return this.startPositionQuery;
                }
            }

            /// <summary>
            /// Get the end position of the database sequence with respect to the database parent
            /// </summary>
            public int EndDatabase
            {
                get
                {
                    return endPositionDatabase;
                }
            }


            /// <summary>
            /// Get the end position of the query sequence with respect to the query parent
            /// </summary>
            public int EndQuery
            {
                get
                {
                    return this.endPositionQuery;
                }
            }

            /// <summary>
            /// Get the total number of gaps. Actually is the number of occurences of the char '-' within either the database or the query sequence.
            /// </summary>
            public int Gaps
            {
                get
                {
                    if (this.gaps == null)
                    {
                        this.gaps = 0;
                        for (int i = 0; i < alignments.Count; i++)
                        {
                            this.gaps += alignments[i].Gaps;
                        }
                    }
                    return gaps.Value;

                }
            }

            /// <summary>
            /// Get the total length of the alignment when gaps are considered
            /// </summary>
            public int LengthAlignmentWithGaps
            {
                get
                {
                    if (lengthWithGaps == null)
                    {
                        this.lengthWithGaps = 0;
                        for (int i = 0; i < alignments.Count; i++)
                        {
                            this.lengthWithGaps += alignments[i].LengthAlignmentWithGaps;
                        }
                    }
                    return lengthWithGaps.Value;
                }
            }


            /// <summary>
            /// Get the total length of the alignment, not considering any gaps
            /// </summary>
            public int LengthAlignmentWithoutGaps
            {
                get
                {
                    if (lengthWithoutGaps == null)
                    {
                        this.lengthWithoutGaps = 0;
                        for (int i = 0; i < alignments.Count; i++)
                        {
                            this.lengthWithoutGaps += alignments[i].LengthAlignmentWithoutGaps;
                        }
                    }
                    return lengthWithoutGaps.Value;
                }
            }

            /// <summary>
            /// Get the total score of the alignment
            /// </summary>
            public float? Score
            {
                get
                {
                    if (this.score == null)
                    {
                        this.score = 0;
                        for (int i = 0; i < alignments.Count; i++)
                        {
                            this.score += alignments[i].Score;
                        }
                    }
                    return score;

                }
                set
                {
                    throw new InvalidOperationException("It is not possible to set the score in the composite mode. The score is automatically calculated as sum of the scores of the subalignments");
                }
            }

            /// <summary>
            /// Get the similsarity of the alignment in percent. Gaps are not considered
            /// </summary>
            public float SimilarityWithoutGaps
            {
                get
                {
                    if (this.similaritywithoutgaps == null)
                    {
                        this.similaritywithoutgaps = 0.0F;
                        float sum = 0.0F;
                        for (int i = 0; i < alignments.Count; i++)
                        {
                            sum += alignments[i].SimilarityWithoutGaps * alignments[i].LengthAlignmentWithoutGaps;
                        }
                        this.similaritywithoutgaps = sum / this.LengthAlignmentWithoutGaps;
                    }
                    return similaritywithoutgaps.Value;


                }
            }

            /// <summary>
            /// Get the similarity of the alignment in percent. Gaps are considered
            /// </summary>
            public float SimilarityWithGaps
            {
                get
                {
                    if (this.similaritywithgaps == null)
                    {
                        this.similaritywithgaps = 0.0F;
                        float sum = 0.0F;
                        for (int i = 0; i < alignments.Count; i++)
                        {
                            sum += alignments[i].SimilarityWithGaps * alignments[i].LengthAlignmentWithGaps;
                        }
                        this.similaritywithgaps = sum / this.LengthAlignmentWithGaps;
                    }
                    return similaritywithgaps.Value;
                }
            }

            /// <summary>
            /// Get the 
            /// </summary>
            /// <param name="index"></param>
            /// <returns></returns>
            public IPairwiseAlignmentContainer this[int index]
            {
                get
                {
                    return this.alignments[index];
                }
            }

            /// <summary>
            /// Get the number of subalignments
            /// </summary>
            public int Count_SubAlignments
            {
                get
                {
                    return alignments.Count;
                }
            }

            /// <summary>
            /// Gets or sets the instance of the ISignificator which will be used to determine whether or not this alignment is significant
            /// </summary>
            public Bio.Blast.ISignificator Significator
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

            /// <summary>
            /// Is this alignment significant.
            /// The given ISignificator is used (strategy pattern)
            /// </summary>
            /// <returns></returns>
            public bool IsSignificant()
            {

                return this.significator.IsSignificant(this);
            }


            /// <summary>
            /// Return the length of this alignment. Is length of the alignment with respect to the database sequence without gaps
            /// </summary>
            public long Length
            {
                get
                {
                    if (this.length == null)
                    {
                        this.length = 0;
                        for (int i = 0; i < alignments.Count; i++)
                        {
                            this.length += alignments[i].Length;
                        }
                    }
                    return this.length.Value;

                }
            }




            /// <summary>
            /// Get the start position of this alignment.
            /// </summary>
            /// <exception cref="InvalidOperationException">Set property is not implemented</exception>
            public int? Start
            {
                get
                {
                    return this.StartDatabase;
                }
                set
                {
                    throw new InvalidOperationException("Start property can not be assigned after instantiation; It is read only");
                }
            }

            /// <summary>
            /// Get the end position of the alignmetn
            /// </summary>
            /// <exception cref="InvalidOperationException">Set property is not implemented</exception>
            public int? End
            {
                get
                {
                    return this.EndDatabase;
                }
                set
                {
                    throw new InvalidOperationException("End property can not be assigned after instantiation; It is read only");
                }
            }

            /// <summary>
            /// Is this alignment the root. Returns false, alignments can never be the root
            /// </summary>
            /// <exception cref="InvalidOperationException">Set property is not implemented</exception>
            public bool? IsRoot
            {
                get
                {
                    return false;
                }
                set
                {
                    throw new InvalidOperationException("An alignment is always part of database sequence and can therefore never be the root");
                }
            }

            /// <summary>
            /// Get or set the name of the database parent
            /// </summary>
            public string ParentName
            {
                get
                {
                    return this.DatabaseParent;
                }
                set
                {
                    this.DatabaseParent = value;
                }
            }





            /// <summary>
            /// Get the total number of transitions
            /// </summary>
            public int Transitions
            {
                get
                {
                    if (this.transitions == null)
                    {
                        transitions = 0;
                        foreach (PairwiseNucleotideSequenceAlignment pw in alignments)
                        {
                            transitions += pw.Transitions;
                        }
                    }
                    return transitions.Value;
                }
            }



            /// <summary>
            /// Get the total number of transversions
            /// </summary>
            public int Transversions
            {
                get
                {
                    if (this.transversions == null)
                    {
                        transversions = 0;
                        foreach (PairwiseNucleotideSequenceAlignment pw in alignments)
                        {
                            transversions += pw.Transversions;
                        }
                    }
                    return transversions.Value;
                }
            }






            /// <summary>
            /// Count the number of long gaps found within both the database and the query sequence;
            /// Only gaps will be considered which exceed the specified minimum gap length (or equal);
            /// </summary>
            /// <param name="minGapLength">minimum length of the gap</param>
            /// <returns></returns>
            public int Count_LongGaps(int minGapLength)
            {
                return this.Count_LongGaps_Database(minGapLength) + this.Count_LongGaps_Query(minGapLength);
            }

            public int Count_LongGaps_Database(int minGapLength)
            {
                int count = 0;
                foreach (PairwiseNucleotideSequenceAlignment pw in alignments)
                {
                    count += pw.Count_LongGaps_Database(minGapLength);
                }
                if (alignments.Count > 1)
                {
                    for (int i = 1; i < alignments.Count; i++)
                    {
                        if (alignments[i].StartQuery - alignments[i - 1].EndQuery - 1 >= minGapLength) count++;
                    }
                }
                return count;
            }

            /// <summary>
            /// Count the number of long gaps found in the query sequence which exceed the specified minimum length of the gap.
            /// Considers both the gaps in the subalignments and the gaps between the subalignments
            /// </summary>
            /// <param name="minGapLength">minimum length of the gap</param>
            /// <returns></returns>
            public int Count_LongGaps_Query(int minGapLength)
            {
                int count = 0;
                foreach (PairwiseNucleotideSequenceAlignment pw in alignments)
                {
                    count += pw.Count_LongGaps_Query(minGapLength);
                }
                if (alignments.Count > 1)
                {
                    for (int i = 1; i < alignments.Count; i++)
                    {
                        if (alignments[i].StartDatabase - alignments[i - 1].EndDatabase - 1 >= minGapLength) count++;
                    }
                }
                return count;
            }


            /// <summary>
            /// Get a collection of gaps
            /// </summary>
            /// <param name="minGapLength"></param>
            /// <returns></returns>
            public List<GapInPairwiseNucleotideSequenceAlignment> LongGaps_Query(int minGapLength)
            {
                List<GapInPairwiseNucleotideSequenceAlignment> gaps = new List<GapInPairwiseNucleotideSequenceAlignment>();

                foreach (PairwiseNucleotideSequenceAlignment pw in alignments)
                {
                    gaps.AddRange(pw.LongGaps_Query(minGapLength));
                }

                if (alignments.Count > 1)
                {
                    for (int i = 1; i < alignments.Count; i++)
                    {
                        if (alignments[i].StartDatabase - alignments[i - 1].EndDatabase - 1 >= minGapLength)
                            gaps.Add(new GapInPairwiseNucleotideSequenceAlignment(minGapLength, alignments[i - 1].EndDatabase, alignments[i].StartDatabase, alignments[i - 1].EndQuery, alignments[i].StartQuery, this));
                    }
                }
                return gaps;
            }


            /// <summary>
            /// 
            /// </summary>
            /// <param name="minGapLength"></param>
            /// <returns></returns>
            public List<GapInPairwiseNucleotideSequenceAlignment> LongGaps_Database(int minGapLength)
            {
                List<GapInPairwiseNucleotideSequenceAlignment> gaps = new List<GapInPairwiseNucleotideSequenceAlignment>();


                foreach (PairwiseNucleotideSequenceAlignment pw in alignments)
                {
                    gaps.AddRange(pw.LongGaps_Database(minGapLength));
                }


                if (alignments.Count > 1)
                {
                    for (int i = 1; i < alignments.Count; i++)
                    {
                        if (alignments[i].StartQuery - alignments[i - 1].EndQuery - 1 >= minGapLength)
                            gaps.Add(new GapInPairwiseNucleotideSequenceAlignment(minGapLength, alignments[i - 1].EndDatabase, alignments[i].StartDatabase, alignments[i - 1].EndQuery, alignments[i].StartQuery, this));

                    }
                }
                return gaps;
            }

            public List<GapInPairwiseNucleotideSequenceAlignment> LongGaps(int minGapLength)
            {
                List<GapInPairwiseNucleotideSequenceAlignment> gaps = new List<GapInPairwiseNucleotideSequenceAlignment>();
                gaps.AddRange(LongGaps_Query(minGapLength));
                gaps.AddRange(LongGaps_Database(minGapLength));
                return gaps;
            }






            public int? LengthDatabaseParent
            {
                get
                {
                    if (alignments.Count > 0) return this.alignments[0].LengthDatabaseParent;
                    else return -1;
                }
                set
                {

                    foreach (PairwiseNucleotideSequenceAlignment pw in alignments)
                    {
                        pw.LengthDatabaseParent = value;
                    }
                }
            }

            public int? LengthQueryParent
            {
                get
                {
                    if (alignments.Count > 0) return this.alignments[0].LengthQueryParent;
                    else return -1;
                }
                set
                {
                    foreach (PairwiseNucleotideSequenceAlignment pw in alignments)
                    {
                        pw.LengthQueryParent = value;
                    }
                }
            }

            public List<IPairwiseAlignmentContainer> SubAlignments
            {
                get
                {
                    return this.alignments.ConvertAll<IPairwiseAlignmentContainer>(new Converter<PairwiseNucleotideSequenceAlignment, IPairwiseAlignmentContainer>
                (delegate(PairwiseNucleotideSequenceAlignment toConvert) { return (IPairwiseAlignmentContainer)toConvert; }));
                }
            }





            public PairwiseNucleotideSequenceAlignment GetAlignmentCoveringDatabasePosition(int databasePosition)
            {
                foreach (PairwiseNucleotideSequenceAlignment pw in this.alignments)
                {
                    if (pw.StartDatabase <= databasePosition && pw.EndDatabase >= databasePosition) return pw;
                }
                return null;
            }


            public PairwiseNucleotideSequenceAlignment GetAlignmentCoveringQueryPosition(int queryPosition)
            {
                foreach (PairwiseNucleotideSequenceAlignment pw in this.alignments)
                {
                    if (pw.StartQuery <= queryPosition && pw.EndQuery >= queryPosition) return pw;
                }
                return null;
            }

            public object Tag
            {
                get
                {
                    return this.tag;
                }
                set
                {
                    this.tag = value;
                }
            }





            public int Count_TotalQueryHits
            {
                get
                {
                    return this.countTotalQueryHits;
                }
                set
                {
                    this.countTotalQueryHits = value;
                }
            }






            public string FeatureName
            {
                get { return "Composite pairwise nucleotide sequence alignment"; }
            }

        }


        /// <summary>
        /// Represents a gap in a pairwise nucleotide sequence alignments,
        /// usually such gaps are discoverd for a specified minimum gap length;
        /// </summary>
        public class GapInPairwiseNucleotideSequenceAlignment : IPositionable
        {
            //mandatory parameters
            private int startPositionDatabase;
            private int minimumGapLength;
            private int startPositionQuery;
            private int endPositionDatabase;
            private int endPositionQuery;
            private IPairwiseNucleotideSequenceAlignmentContainer parentAlignment;

            //calculated paramters
            private int? gapLengthQuery;
            private int? gapLengthDatabase;



            //facultative parameters
            private string databaseParent;
            private string queryParent;



            public GapInPairwiseNucleotideSequenceAlignment(int minimumGapLength, int startPositionDatabase, int endPositionDatabase, int startPositionQuery, int endPositionQuery, IPairwiseNucleotideSequenceAlignmentContainer parentAlignment)
            {
                this.minimumGapLength = minimumGapLength;
                this.startPositionDatabase = startPositionDatabase;
                this.startPositionQuery = startPositionQuery;
                this.endPositionDatabase = endPositionDatabase;
                this.endPositionQuery = endPositionQuery;
                this.parentAlignment = parentAlignment;
            }


            public int StartPositionDatabase
            {
                get
                {
                    return this.startPositionDatabase;
                }
            }
            public int StartPositionQuery
            {
                get
                {
                    return this.startPositionQuery;
                }
            }
            public int EndPositionDatabase
            {
                get
                {
                    return this.endPositionDatabase;
                }
            }
            public int EndPositionQuery
            {
                get
                {
                    return this.endPositionQuery;
                }
            }
            public int GapLengthDatabase
            {
                get
                {
                    if (this.gapLengthDatabase == null)
                    {
                        this.gapLengthDatabase = endPositionQuery - startPositionQuery - 1;
                    }
                    return gapLengthDatabase.Value;

                }
            }
            public int GapLengthQuery
            {
                get
                {
                    if (this.gapLengthQuery == null)
                    {
                        this.gapLengthQuery = endPositionDatabase - startPositionDatabase - 1;
                    }
                    return this.gapLengthQuery.Value;
                }

            }

            public bool QueryGap
            {
                get
                {
                    if (this.GapLengthQuery > this.minimumGapLength) return true;
                    else return false;
                }
            }


            public bool DatabaseGap
            {
                get
                {
                    if (this.GapLengthDatabase >= this.minimumGapLength) return true;
                    else return false;
                }
            }

            public string DatabaseParent
            {
                get
                {
                    return this.databaseParent;
                }
                set
                {
                    this.databaseParent = value;
                }
            }

            public string QueryParent
            {
                get
                {
                    return this.queryParent;
                }
                set
                {
                    this.queryParent = value;
                }
            }





            int? IPositionable.Start
            {
                get
                {
                    return this.StartPositionDatabase;
                }
                set
                {
                    throw new Exception("The method or operation is not implemented.");
                }
            }

            public int? End
            {
                get
                {
                    return this.EndPositionDatabase;
                }
                set
                {
                    throw new Exception("The method or operation is not implemented.");
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
                    throw new Exception("The method or operation is not implemented.");
                }
            }

            public string ParentName
            {
                get
                {
                    return this.DatabaseParent;
                }
                set
                {
                    this.DatabaseParent = value;
                }
            }





            public long Length
            {
                get { throw new Exception("The method or operation is not implemented."); }
            }


        



            public string FeatureName
            {
                get { return "Gap in a nucleotide sequence alignment"; }
            }


        }


        public interface IAlignmentFormater
        {

            ISequenceContainer GetSimilaritySequence(PairwiseAlignment pwAlignment);
        }


        /// <summary>
        /// Formater class for AlignmentOfTwoSequences
        /// Creates and formats a string representations for a given alignment in the way as identified in the Water program of Emboss
        /// Identical characters are represented through a '|' mismatches through a '.' and gaps through ' '
        /// </summary>
    public class PairwiseAlignmentFormater_WaterEmboss : IAlignmentFormater
    {

        Dictionary<char, bool> badChars;


        public PairwiseAlignmentFormater_WaterEmboss()
        {

            badChars = new Dictionary<char, bool>();
            badChars.Add('n', true);
            badChars.Add('N', true);
            badChars.Add('X', true);
            badChars.Add('x', true);
            badChars.Add('-', true);
        }


        public ISequenceContainer GetSimilaritySequence(PairwiseAlignment pwAlignment)
        {

            ISequenceContainer similaritySequence = SequenceFactory.GetDefaultSequence();
            ISequenceContainer data = pwAlignment.DatabaseSequence;
            ISequenceContainer query = pwAlignment.QuerySequence;

            similaritySequence.Capacity = pwAlignment.Length;

            for (int i = 0; i < pwAlignment.Length; i++)
            {
                if (badChars.ContainsKey(query[i]) || badChars.ContainsKey(data[i]))
                {
                    similaritySequence.Append(' ');
                }
                else similaritySequence.Append(data[i] == query[i] ? '|' : '.');

            }
            return similaritySequence;
        }
    }
            


        


        /// <summary>
        /// Formater class for AlignmentOfTwoSequences
        /// Creates and formats a string representations for a given alignment in the way as identified in the Water program of Emboss
        /// Identical characters are represented through a '|' mismatches through a '.' and gaps through ' '
        /// </summary>
    public class PairwiseAlignmentFormater_Blast : IAlignmentFormater
    {

        Dictionary<char, bool> badChars;

        public PairwiseAlignmentFormater_Blast()
        {

            badChars = new Dictionary<char, bool>();
            badChars.Add('n', true);
            badChars.Add('N', true);
            badChars.Add('X', true);
            badChars.Add('x', true);
            badChars.Add('-', true);
        }




        public ISequenceContainer GetSimilaritySequence(PairwiseAlignment pwAlignment)
        {

            ISequenceContainer similaritySequence = SequenceFactory.GetDefaultSequence();
            similaritySequence.Capacity = pwAlignment.Length;
            ISequenceContainer data = pwAlignment.DatabaseSequence;
            ISequenceContainer query = pwAlignment.QuerySequence;
            for (int i = 0; i < pwAlignment.Length; i++)
            {
                if (badChars.ContainsKey(query[i]) || badChars.ContainsKey(data[i]))
                {
                    similaritySequence.Append(' ');
                }
                else similaritySequence.Append(data[i] == query[i] ? '|' : ' ');

            }
            return similaritySequence;
        }
    }

        





        /// <summary>
        /// Basic container for the alignment of two sequences, irrespective if nucleotide or protein sequence
        /// </summary>
        public class PairwiseAlignment
        {
            private ISequenceContainer databaseSequence = null;
            private ISequenceContainer querySequence = null;
            private ISequenceContainer similaritySequence = null;


            private IAlignmentFormater formater = new PairwiseAlignmentFormater_WaterEmboss();





            public PairwiseAlignment(ISequenceContainer databaseSequence, ISequenceContainer querySequence)
            {
                if (databaseSequence.Length != querySequence.Length) throw new InvalidOperationException("Unable to initialize TwoSequenceAlignment, the two sequences have to be of equal length");
                this.databaseSequence = databaseSequence;
                this.querySequence = querySequence;
            }





            /// <summary>
            /// Returns the alignment for a given index as byte array
            /// Examples {'A',' ','T'} or {'A','|','A'}
            /// position[0]= database character
            /// position[1]= character indicating equality
            /// position[2]= query character
            /// </summary>
            public char[] this[int index]
            {

                get
                {
                    char[] c = new char[3];
                    c[0] = DatabaseSequence[index];
                    c[1] = SimilaritySequence[index];
                    c[2] = QuerySequence[index];
                    return c;
                }
            }


            /// <summary>
            /// Returns the ISequence of the DatabaseSequence with the character '-' marking insertions in the QuerySequence with
            /// respect to the DatabaseSequence
            /// eg. AATGCTAG--TGGCTTC
            /// </summary>
            public ISequenceContainer DatabaseSequence
            {
                get
                {
                    return databaseSequence;
                }
            }

            public IAlignmentFormater AlignmentFormater
            {
                get
                {
                    return this.formater;
                }
                set
                {
                    this.formater = value;
                    this.similaritySequence = null;
                }
            }


            /// <summary>
            /// Returns the ISequence of the DatabaseSequence without the characters '-', which depict indels
            /// </summary>
            public ISequenceContainer DatabaseSequence_WithoutIndels
            {
                get
                {
                    ISequenceContainer ret = SequenceFactory.GetDefaultSequence();
                    for (int i = 0; i < DatabaseSequence.Length; i++)
                    {
                        if (DatabaseSequence[i] != '-') ret.Append(DatabaseSequence[i]);
                    }
                    return ret;
                }
            }
            /// <summary>
            /// Returns the sequence depicting similarity between the DatabaseSequence and the QuerySequence
            /// eg. ||||  |||  ||||||
            /// </summary>
            public ISequenceContainer SimilaritySequence
            {
                get
                {
                    if (similaritySequence == null) similaritySequence = formater.GetSimilaritySequence(this);
                    return similaritySequence;
                }
            }

            /// <summary>
            /// Returns the ISequence of a QuerySequence, with the characters '-' marking deletions in the Query sequence 
            /// with respect to the Database Sequence
            /// eg. AATCTC-ACTC-TAG
            /// </summary>
            public ISequenceContainer QuerySequence
            {
                get
                {
                    return querySequence;
                }
            }


            /// <summary>
            /// Returns the ISequence of a QuerySequence, without the characters '-' which are usually used for depicting indels
            /// eg. AATGTCTCAATGTTTG
            /// </summary>
            public ISequenceContainer QuerySequence_WithoutIndels
            {
                get
                {
                    ISequenceContainer ret = SequenceFactory.GetDefaultSequence();
                    for (int i = 0; i < QuerySequence.Length; i++)
                    {
                        if (QuerySequence[i] != '-') ret.Append(QuerySequence[i]);
                    }
                    return ret;
                }
            }


            /// <summary>
            /// Returns a sub-alignment
            /// </summary>
            public PairwiseAlignment GetSubAlignment(int startPosition, int length)
            {
                return new PairwiseAlignment(databaseSequence.SubSequence(startPosition, length), querySequence.SubSequence(startPosition, length));

            }

            public int Length
            {
                get
                {
                    return databaseSequence.Length;
                }
            }

            /// <summary>
            /// Returns the alignment of two sequences in the form specified by the formater
            /// </summary>
            /// <returns></returns>
            public override string ToString()
            {
                StringBuilder sb = new StringBuilder();
                sb.Capacity = (databaseSequence.Length + 3) * 3;
                sb.Append(DatabaseSequence + Environment.NewLine);
                sb.Append(SimilaritySequence + Environment.NewLine);
                sb.Append(QuerySequence + Environment.NewLine);
                return sb.ToString();
            }


        }

    public static class PairwiseAlignmentUtility
    {
        public static PairwiseAlignment ReverseAlignment(PairwiseAlignment alignment)
        {
            if (alignment == null) return null;
            ISequenceContainer data = alignment.DatabaseSequence;
            ISequenceContainer query = alignment.QuerySequence;
            //Reverse the individual sequences
            return new PairwiseAlignment(SequenceUtility.ReverseSequence(data), SequenceUtility.ReverseSequence(query));
        }
    }







}
