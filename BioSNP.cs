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
using Bio.Seq;
using Bio.Alignment;

namespace Bio.SNP
{

    public class RefSNP : IPositionable
    {
        //532	2	0	1	1	1	14	NW_925539	1	HsCraAADB02_477	19457932	19515810	PNN	0.19848	0.244634	 	15	1	1	36	126	Celera
        //12      Position of RefSNP in chromosome coordinates (used to order report)
        // Locations are specified in NCBI sequence location convention where:
        //      x, a single number, indicates a feature at base position x
        //      x..y, denotes a feature that spans from x to y inclusive.
        //      x^y, denotes a feature that is inserted between bases x and y
        private int refSNPid;
        private int mapweight;
        private int snp_type;
        private int chromosomeHits;
        private int contigHits;
        private int genomeHits;
        private string contigAccession;
        private int contigVersionNumber;
        private string contigID;
        private string chromosome;
        private int? start;
        private int? start_contig;
        private string parentName;
        private string[] genes;
        private double? averageHeterozygosity;
        private double? stdDevAvHet;
        private double? preal;
        private int validationStatus;
        private bool genotypeAvailable;
        private bool linkOutAvailable;
        private int dbsnpIdCreation;
        private int dbSNPmostrecentchange;
        private string mappedToReference;

        public RefSNP(int refSNPid, int mapweight, int snpType, int chromosomeHits, int contigHits, int genomeHits, string chromosome,
            string contigAccession, int contigVersionNumber, string contigID, int? startposition_contig, int? startPosition,
            string[] genes, double? averageHeterozygosity, double? stdDev_averageHeterozygosity, double? p_realSNP, int validationStatus,
            bool genotypeAvailable, bool linkOutAvailable, int dbSNPIdCreation, int dbSNPmostRecentChange, string mappedToReference)
        {
            this.refSNPid = refSNPid;
            this.mapweight = mapweight;
            this.snp_type = snpType;
            this.chromosomeHits = chromosomeHits;
            this.contigHits = contigHits;
            this.genomeHits = genomeHits;
            this.contigAccession = contigAccession;
            this.contigVersionNumber = contigVersionNumber;
            this.contigID = contigID;
            this.chromosome = chromosome;
            this.start = startPosition;
            this.start_contig = startposition_contig;
            this.genes = genes;
            this.averageHeterozygosity = averageHeterozygosity;
            this.stdDevAvHet = stdDev_averageHeterozygosity;
            this.preal = p_realSNP;
            this.validationStatus = validationStatus;
            this.genotypeAvailable = genotypeAvailable;
            this.linkOutAvailable = linkOutAvailable;
            this.dbsnpIdCreation = dbSNPIdCreation;
            this.dbSNPmostrecentchange = dbSNPmostRecentChange;
            this.mappedToReference = mappedToReference;
        }

        public RefSNP()
        {
        }

        #region Properties
        /// <summary>
        /// ID of the SNP
        /// </summary>
        public int RefSNPid
        {
            get { return refSNPid; }
            set { refSNPid = value; }
        }
        /// <summary>
        /// 1 = Unmapped
        /// 2 = Mapped to single position in genome
        /// 3 = Mapped to 2 positions on a single chromosome
        /// 4 = Mapped to 3-10 positions in genome (possible paralog hits)
        /// 5 = Mapped to >10 positions in genome
        /// </summary>
        public int Mapweight
        {
            get { return mapweight; }
            set { mapweight = value; }
        }

        /// <summary>
        /// 0 = Not withdrawn.
        /// 1 = Withdrawn. There are several reasons for withdrawn, the
        /// withdrawn status is fully defined in the asn1, flatfile,
        /// and XML descriptions of the RefSNP. See /specs/docsum_2005.asn
        /// for a full definition of snp-type values.
        /// </summary>
        public int SNP_type
        {
            get { return snp_type; }
            set { snp_type = value; }
        }

        /// <summary>
        /// Total number of chromosomes hit by this RefSNP during mapping
        /// </summary>
        public int ChromosomeHits
        {
            get { return chromosomeHits; }
            set { chromosomeHits = value; }
        }

        /// <summary>
        /// Total number of contigs hit by this RefSNP during mapping
        /// </summary>
        public int ContigHits
        {
            get { return contigHits; }
            set { contigHits = value; }
        }


        /// <summary>
        /// Total number of hits to genome by this RefSNP during mapping          
        /// </summary>
        public int GenomeHits
        {
            get { return genomeHits; }
            set { genomeHits = value; }
        }


        /// <summary>
        /// Chromosome number
        /// </summary>
        public string Chromosome
        {
            get { return chromosome; }
            set { chromosome = value; }
        }

        /// <summary>
        /// Contig accession for this hit to genome
        /// </summary>
        public string ContigAccesssion
        {
            get { return contigAccession; }
            set { contigAccession = value; }
        }



        /// <summary>
        /// Version number of contig accession for this hit to genome
        /// </summary>
        public int ContigVersionNumber
        {
            get { return contigVersionNumber; }
            set { contigVersionNumber = value; }
        }



        /// <summary>
        /// Contig ID for this hit to genome
        /// </summary>
        public string ContigID
        {
            get { return contigID; }
            set { contigID = value; }
        }



        /// <summary>
        /// Position of RefSNP in contig coordinates
        /// </summary>
        public int? Start_Contig
        {
            get { return start_contig; }
            set { start_contig = value; }
        }

        /// <summary>
        ///  12      Position of RefSNP in chromosome coordinates (used to order report)
        ///Locations are specified in NCBI sequence location convention where:
        /// </summary>
        public int? Start
        {
            get
            {
                return this.start;
            }
            set
            {
                if (value != null)
                {
                    this.start = value.Value;
                }
            }
        }

        public int? End
        {
            get
            {
                return this.start;
            }
            set
            {
                throw new InvalidOperationException("The end position of a SNP is the same as the start position; setting of this property is not allowed");
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
                throw new InvalidOperationException("A SNP can never be the root of a sequence, this property will therefore always return false! Setting of this property is not valid");
            }
        }

        public string ParentName
        {
            get
            {
                if (parentName != null && parentName != "") return parentName;
                else if (chromosome != null && chromosome != "") return chromosome;
                else return null;
            }
            set
            {
                this.parentName = value;
            }
        }
        public long Length
        {
            get
            {
                return 1;
            }
            set
            {
                throw new InvalidOperationException("Length of a SNP is always exactly: 1; setting this property to a value is not allowed");
            }
        }

        /// <summary>
        /// Genes at this same position on the chromosome
        /// </summary>
        public string[] Genes
        {
            get
            {
                return genes;
            }
            set
            {
                this.genes = value;
            }
        }


        /// <summary>
        /// Average heterozygosity of this RefSNP
        /// </summary>
        public double? AverageHeterozygosity
        {
            get { return averageHeterozygosity; }
            set { averageHeterozygosity = value; }
        }

        /// <summary>
        /// Standard error of average heterozygosity
        /// </summary>
        public double? StdDev_AverageHeterozygosity
        {
            get { return stdDevAvHet; }
            set { stdDevAvHet = value; }
        }


        /// <summary>
        /// Maximum reported probability that RefSNP is real. (For computationally-
        ///predicted submissions)
        /// </summary>
        public double? P_real
        {
            get { return preal; }
            set { preal = value; }
        }

        /// <summary>
        /// Validated status
        /// 0 = No validation information
        /// 1 = Cluster has 2+ submissions, with 1+ submission assayed 
        ///     with a non-computational method
        /// 2 = At least one subsnp in cluster has frequency data submitted
        /// 3 = Non-computational method in cluster and frequency data present
        /// 4 = At lease one subsnp in cluster has been experimentally 
        ///     validated by submitter
        ///     for other validation status values, please see:
        ///     a href="ftp://ftp.ncbi.nih.gov/snp/database/organism_shared_data/SnpValidationCode.bcp.gz">ftp://ftp.ncbi.nih.gov/snp/database/organism_shared_data
        ///     /SnpValidationCode.bcp.gz
        /// </summary>
        public int ValidationStatus
        {
            get { return validationStatus; }
            set { validationStatus = value; }
        }


        /// <summary>
        /// Genotypes available in dbSNP for this RefSNP
        /// </summary>
        public bool GenotypeAvailable
        {
            get { return genotypeAvailable; }
            set { genotypeAvailable = value; }
        }

        /// <summary>
        ///   19      Linkout available to submitter website for further data on the RefSNP     
        /// </summary>
        public bool LinkoutAvailable
        {
            get { return linkOutAvailable; }
            set { linkOutAvailable = value; }
        }


        /// <summary>
        /// dbSNP build ID when the refSNP was first created (i.e. the create date)
        /// </summary>
        public int DbSNPID_FirstCreated
        {
            get { return dbsnpIdCreation; }
            set { dbsnpIdCreation = value; }
        }

        /// <summary>
        /// dbSNP build ID of the most recent change to the refSNP cluster. The
        /// date of the change is represented by the build ID which has an
        ////approximate date/time associated with it. (see:
        /// http://www.ncbi.nlm.nih.gov/projects/SNP/buildhistory.cgi)
        /// </summary>
        public int DbSNPID_MostRecentChange
        {
            get { return dbSNPmostrecentchange; }
            set { dbSNPmostrecentchange = value; }
        }


        /// <summary>
        /// Mapped to a reference or alternate (e.g. Celera) assembly 
        /// </summary>
        public string MappedToReference
        {
            get { return mappedToReference; }
            set { mappedToReference = value; }
        }



















        public string FeatureName
        {
            get {return "RefSNP"; }
        }

        #endregion





    }

    public class SNPUtility
    {
        /// <summary>
        /// Extract SNPLoci from a given collection of SNPs. 
        /// </summary>
        /// <param name="snps">the collection of SNPs which will be used to obtain the collection of SNP sites</param>
        /// <returns></returns>
        public static List<SNPSite> GetSNPSites(List<SNP> snps)
        {
            List<SNPSite> sites = new List<SNPSite>();

            List<SNP> temp = new List<SNP>();
            snps.Sort(new Sort.SortSNPs_DatabaseName_DatabasePosition());
            if (snps.Count > 0)
            {
                temp.Add(snps[0]);
            }
            else return sites; //if the collection, snps, does not contain an entry return an empty collection.


            for (int i = 1; i < snps.Count; i++)
            {
                SNP lastSNP = temp[temp.Count - 1];

                if (lastSNP.DatabaseParent == snps[i].DatabaseParent && lastSNP.SiteBasis_Database == snps[i].SiteBasis_Database && lastSNP.SiteIndelShift_Database == snps[i].SiteIndelShift_Database)
                {
                    //Increase Count
                    temp.Add(snps[i]);
                }
                else
                {

                    sites.Add(new SNPSite(temp));
                    temp = new List<SNP>();
                    temp.Add(snps[i]);
                }

            }
            sites.Add(new SNPSite(temp)); //Also add the last snps

            /*
                            //Also set the total count
                            foreach (SNPSite sn in sites)
                            {
                                sn.Number_SNPsAtDatabaseParent = snps.Count;
                            }
             */

            return sites; //Return the results


        }


        /// <summary>
        /// Get a subset of a given collection of SNPs.
        /// Only SNPs being located at one of the specified locations are returned
        /// Ingenious solution
        /// </summary>
        /// <param name="snps">the collection of snps which will be the basis for creating the subset</param>
        /// <param name="snpSites">the collection of SNPsites, the SNPs must have a position at one of the given sites</param>
        /// <returns></returns>
        public static List<SNP> GetSNPsHavingSNPSite(List<SNP> snps, List<SNPSite> snpSite)
        {
            List<SNP> ret = new List<SNP>();
            if (snps.Count < 1) return ret;


            //DatabaseParent - SiteBasis - SiteIndelShift
            Dictionary<string, Dictionary<int, Dictionary<int, byte>>> dict = new Dictionary<string, Dictionary<int, Dictionary<int, byte>>>();
            foreach (SNPSite locus in snpSite) //Very cool create a Dictionary - Dictionary 
            {
                //The first key is the database parent
                if (dict.ContainsKey(locus.DatabaseParent))
                {
                    //The second key is the database position
                    Dictionary<int, Dictionary<int, byte>> tempDict = dict[locus.DatabaseParent];

                    if (tempDict.ContainsKey(locus.SiteBase_Database))
                    {
                        Dictionary<int, byte> tDict = tempDict[locus.SiteBase_Database];

                        if (!tDict.ContainsKey(locus.SiteIndelShift_Database))
                        {
                            //Database parent exists - SNP site basis exists - the indel shift not
                            tDict.Add(locus.SiteIndelShift_Database, (byte)'0');
                        }

                    }
                    else //Parent exits - SNP site not
                    {
                        // Create a new SNP site for an already existing parent
                        Dictionary<int, byte> temp = new Dictionary<int, byte>();
                        temp.Add(locus.SiteIndelShift_Database, (byte)'0');
                        tempDict.Add(locus.SiteBase_Database, temp);

                    }
                }
                else //Not even the parent exists
                {
                    //Create a totaly new entry
                    Dictionary<int, Dictionary<int, byte>> tempDict = new Dictionary<int, Dictionary<int, byte>>();

                    ///First the indelshift
                    Dictionary<int, byte> temp2 = new Dictionary<int, byte>();
                    temp2.Add(locus.SiteIndelShift_Database, (byte)'0');

                    //Then the site position
                    tempDict.Add(locus.SiteBase_Database, temp2);

                    //Finally the database parent
                    dict.Add(locus.DatabaseParent, tempDict);
                }
            }


            ///Check each snp, if having the dictionary contains the same databaseParent, than the same SNPSite basis and than the same SNPSiteIndelShift
            foreach (SNP snp in snps)
            {

                if (dict.ContainsKey(snp.DatabaseParent) && dict[snp.DatabaseParent].ContainsKey(snp.SiteBasis_Database) && dict[snp.DatabaseParent][snp.SiteBasis_Database].ContainsKey(snp.SiteIndelShift_Database)) ret.Add(snp);
            }

            return ret;
        }


        /*
         *             /// <summary>
        /// Get a subset of a given collection of SNPs.
        /// Only SNPs being located at one of the specified locations are returned
        /// </summary>
        /// <param name="snps">the collection of snps which will be the basis for creating the subset</param>
        /// <param name="snpSites">the collection of SNPsites, the SNPs must have a position at one of the given sites</param>
        /// <returns></returns>
        public static List<SNP> GetSNPsHavingSNPSite(List<SNP> snps, List<SNPSite> snpSite)
        {
            List<SNP> ret = new List<SNP>();
            if (snps.Count < 1) return ret;

            Dictionary<string, Dictionary<string, byte>> dict = new Dictionary<string, Dictionary<string, byte>>();
            foreach (SNPSite locus in snpSite) //Very cool create a Dictionary - Dictionary 
            {
                //The first key is the database parent
                if (dict.ContainsKey(locus.DatabaseParent))
                {
                    //The second key is the database position
                    Dictionary<string, byte> tempDict = dict[locus.DatabaseParent];
                    if (!tempDict.ContainsKey(locus.SiteID_Database))
                    {
                        tempDict.Add(locus.SiteID_Database, (byte)'0');
                    }
                }
                else
                {
                    //Create a totaly new entry
                    Dictionary<string, byte> tempDict = new Dictionary<string, byte>();
                    tempDict.Add(locus.SiteID_Database, (byte)'0');
                    dict.Add(locus.DatabaseParent, tempDict);
                }
            }


            ///Check each snp, if having the dictionary contains the same databaseParent and subsequently the same databasePosition
            foreach (SNP snp in snps)
            {

                if (dict.ContainsKey(snp.DatabaseParent) && dict[snp.DatabaseParent].ContainsKey(snp.SiteID_Database)) ret.Add(snp);
            }

            return ret;
        }

        */


        public static Dictionary<string, List<SNP>> GetSNPDictionary_DatabaseName(List<SNP> snps)
        {
            Dictionary<string, List<SNP>> dict = new Dictionary<string, List<SNP>>();
            foreach (SNP s in snps)
            {
                if (dict.ContainsKey(s.DatabaseParent))
                {
                    dict[s.DatabaseParent].Add(s);
                }
                else
                {
                    List<SNP> sl = new List<SNP>();
                    sl.Add(s);
                    dict.Add(s.DatabaseParent, sl);
                }
            }

#if DEBUG
            Dictionary<string, bool> check = new Dictionary<string, bool>();
            foreach (string s in dict.Keys)
            {
                if (check.ContainsKey(s)) throw new Exception();
                else check.Add(s, true);
            }
#endif

            return dict;
        }


        public static Dictionary<string, List<SNP>> GetSNPDictionary_QueryName(List<SNP> snps)
        {
            Dictionary<string, List<SNP>> dict = new Dictionary<string, List<SNP>>();
            foreach (SNP s in snps)
            {
                if (dict.ContainsKey(s.QueryParent))
                {
                    dict[s.QueryParent].Add(s);
                }
                else
                {
                    List<SNP> sl = new List<SNP>();
                    sl.Add(s);
                    dict.Add(s.QueryParent, sl);
                }
            }

#if DEBUG
            Dictionary<string, bool> check = new Dictionary<string, bool>();
            foreach (string s in dict.Keys)
            {
                if (check.ContainsKey(s)) throw new Exception();
                else check.Add(s, true);
            }
#endif

            return dict;
        }


        public static Dictionary<string, List<SNPSite>> GetSNPSiteDictionary_DatabaseName(List<SNP> snps)
        {
            List<SNPSite> loci = SNPUtility.GetSNPSites(snps);
            return SNPSiteUtility.GetSNPSiteDictionary_DatabaseName(loci);
        }

        /// <summary>
        /// Retrieve a collection of SNP sibblings for a given collection of SNPs.
        /// SNP sibblings are SNPs which share the same common parent, eg. gene or chromosome etc
        /// </summary>
        /// <param name="snps"></param>
        /// <returns></returns>
        public static List<SibblingSNPs> GetSNPSibblings(List<SNP> snps)
        {
            Dictionary<string, List<SNPSite>> dict = GetSNPSiteDictionary_DatabaseName(snps);

            List<SibblingSNPs> siblings = new List<SibblingSNPs>();
            foreach (List<SNPSite> loci in dict.Values)
            {
                siblings.Add(new SibblingSNPs(loci));
            }
            return siblings;

        }

    }



    public struct SNPBenchmarks
    {

        private int wrongPositives; //to much

        public int WrongPositives
        {
            get { return wrongPositives; }
            set { wrongPositives = value; }
        }

        private float sensitivity; //TP/(TP+FN)

        public float Sensitivity
        {
            get { return sensitivity; }
            set { sensitivity = value; }
        }

        private float specificity; // TN/(FP+TN)

        public float Specificity
        {
            get { return specificity; }
            set { specificity = value; }
        }

        private int falseNegatives; //missed

        public int FalseNegatives
        {
            get { return falseNegatives; }
            set { falseNegatives = value; }
        }

        private int countDataSubset;

        public int CountDataSubset
        {
            get { return countDataSubset; }
            set { countDataSubset = value; }
        }
        private int countNotdataSubset;

        public int CountNotdataSubset
        {
            get { return countNotdataSubset; }
            set { countNotdataSubset = value; }
        }
        private int countDataAll;

        public int CountDataAll
        {
            get { return countDataAll; }
            set { countDataAll = value; }
        }
        private int countNotdataAll;

        public int CountNotdataAll
        {
            get { return countNotdataAll; }
            set { countNotdataAll = value; }
        }

        private int tn;

        public int Tn
        {
            get { return tn; }
            set { tn = value; }
        }
        private int fn;

        public int Fn
        {
            get { return fn; }
            set { fn = value; }
        }
        private int tp;

        public int Tp
        {
            get { return tp; }
            set { tp = value; }
        }
        private int fp;

        public int Fp
        {
            get { return fp; }
            set { fp = value; }
        }

    }


    /// <summary>
    /// Utility class to handle SNP loci
    /// </summary>
    public class SNPSiteUtility
    {
        public static SNPBenchmarks GetSNPBenchmarks(List<SNPSite> all, List<SNPSite> subset)
        {
            //Sensitivity:  TP / (TP+FN)
            //Specificity:  TN / (FP+TN)
            SNPBenchmarks benchmarks = new SNPBenchmarks();
            int countValidAll = 0;
            int countValidSubset = 0;
            int countAdd = 0;
            foreach (SNPSite s in all)
            {
                if (s.Alleles.ContainsKey(s.DatabaseCharacter))
                {
                    countValidAll++;
                    int data = 0;
                    int notData = 0;
                    foreach (KeyValuePair<char, List<SNP>> kv in s.Alleles)
                    {
                        if (kv.Key == s.DatabaseCharacter)
                        {
                            data = kv.Value.Count;
                        }
                        else if (kv.Value.Count > notData)
                        {
                            notData = kv.Value.Count;
                        }

                    }

                    if (data > notData) benchmarks.CountDataAll++;
                    else if (data == notData) countAdd++;
                    else if (notData > data) benchmarks.CountNotdataAll++;
                    else throw new Exception("Impossible");
                }

            }
            benchmarks.CountDataAll += countAdd / 2;
            benchmarks.CountNotdataAll += countAdd / 2;

            countAdd = 0;
            foreach (SNPSite s in subset)
            {
                if (s.Alleles.ContainsKey(s.DatabaseCharacter))
                {
                    countValidSubset++;
                    int data = 0;
                    int notData = 0;
                    foreach (KeyValuePair<char, List<SNP>> kv in s.Alleles)
                    {
                        if (kv.Key == s.DatabaseCharacter)
                        {
                            data = kv.Value.Count;
                        }
                        else if (kv.Value.Count > notData)
                        {
                            notData = kv.Value.Count;
                        }

                    }


                    if (data > notData) benchmarks.CountDataSubset++;
                    else if (data == notData) countAdd++;
                    else if (notData > data) benchmarks.CountNotdataSubset++;
                    else throw new Exception("Impossible");
                }
            }
            benchmarks.CountDataSubset += countAdd / 2;
            benchmarks.CountNotdataSubset += countAdd / 2;




            benchmarks.Tp = benchmarks.CountNotdataSubset * 2;
            benchmarks.Fp = countValidSubset - benchmarks.Tp;
            if (benchmarks.Fp < 0) benchmarks.Fp = 0;



            benchmarks.Fn = benchmarks.CountNotdataAll * 2 - benchmarks.Tp; //Missed real ones
            benchmarks.Tn = countValidAll - benchmarks.Fn - benchmarks.Tp - benchmarks.Tn;
            if (benchmarks.Tn < 0) throw new Exception("should not happen");

            //Sensitivity:  TP / (TP+FN)
            //Specificity:  TN / (FP+TN)
            benchmarks.Sensitivity = 100.0F * ((float)benchmarks.Tp) / (float)(benchmarks.Tp + benchmarks.Fn);
            benchmarks.Specificity = 100.0F * ((float)benchmarks.Tn) / (float)(benchmarks.Fp + benchmarks.Tn);
            return benchmarks;
        }


        /// <summary>
        /// Retrieve a SNPLocus dictionary in which all SNPLoci are grouped according to their gene ID (database parent),
        /// Such a dictionary for example, is required for the SNP Allele extraction.
        /// </summary>
        /// <param name="snpSites"></param>
        /// <returns></returns>
        public static Dictionary<string, List<SNPSite>> GetSNPSiteDictionary_DatabaseName(List<SNPSite> snpSites)
        {
            Dictionary<string, List<SNPSite>> dict = new Dictionary<string, List<SNPSite>>();

            for (int i = 0; i < snpSites.Count; i++)
            {
                if (dict.ContainsKey(snpSites[i].DatabaseParent))
                {
                    dict[snpSites[i].DatabaseParent].Add(snpSites[i]);
                }
                else
                {
                    List<SNPSite> site = new List<SNPSite>();
                    site.Add(snpSites[i]);
                    dict.Add(snpSites[i].DatabaseParent, site);
                }

            }
            return dict;
        }



        /// <summary>
        /// Retrieve a collection of SNP sibblings for a given collection of SNPSites.
        /// SNP sibblings are SNPs which share the same parent, eg. gene or chromosome etc
        /// </summary>
        /// <param name="snps"></param>
        /// <returns></returns>
        public static List<SibblingSNPs> GetSNPSibblings(List<SNPSite> sites)
        {
            Dictionary<string, List<SNPSite>> dict = SNPSiteUtility.GetSNPSiteDictionary_DatabaseName(sites);

            List<SibblingSNPs> siblings = new List<SibblingSNPs>();
            foreach (List<SNPSite> loci in dict.Values)
            {
                siblings.Add(new SibblingSNPs(loci));
            }
            return siblings;

        }


        /// <summary>
        /// Retrieve all SNPs mapping to a particular SNP site
        /// </summary>
        /// <param name="loci"></param>
        /// <returns></returns>
        public static List<SNP> GetSNPs_Valid(List<SNPSite> site)
        {
            List<SNP> ret = new List<SNP>();
            foreach (SNPSite locus in site)
            {
                ret.AddRange(locus.SNPs_Valid);
            }
            return ret;
        }

        /// <summary>
        /// Retrieve all SNPs mapping to a particular SNP site
        /// </summary>
        /// <param name="loci"></param>
        /// <returns></returns>
        public static List<SNP> GetSNPs_ValidAndInvalid(List<SNPSite> site)
        {
            List<SNP> ret = new List<SNP>();
            foreach (SNPSite locus in site)
            {
                ret.AddRange(locus.SNPs_Valid);
                ret.AddRange(locus.SNPs_Invalid);
            }
            return ret;
        }





    }

    /// <summary>
    /// Represents a SNP site, i.e. a position in the database sequence where a SNP has been discovered;
    /// This class additionally represents a collection of all SNPs at this site
    /// </summary>
    public class SNPSite : IPositionable
    {
        public SNPSite(string databaseParent, string siteID_database)
        {
            this.databaseParent = databaseParent;
            this.SetSiteIDDatabase(siteID_database);
            this.all = new List<SNP>();
            this.invalid = new List<SNP>();
            this.allAlleles = new Dictionary<char, List<SNP>>();
        }
        public SNPSite(List<SNP> snps)
        {
            //if (snps.Count < 1) throw new InvalidDataException("It is not valid to create an empty SNP site; Number of the SNPs must be larger than zero");
            this.all = new List<SNP>();
            this.invalid = new List<SNP>();

            foreach (SNP s in snps)
            {
                if (s.Valid) all.Add(s);
                else invalid.Add(s);
            }

            ///Instantiate the member variables
            if (snps.Count > 0)
            {
                this.databaseParent = snps[0].DatabaseParent;
                this.siteBase_Database = snps[0].SiteBasis_Database;
                this.siteIndelShift_Database = snps[0].SiteIndelShift_Database;

                this.snpAlleles = snps[0].IsFullSNP;
                this.databaseCharacter = snps[0].DatabaseChar;
            }
            else
            {
                databaseParent = null;
                snpAlleles = false;
                this.siteIndelShift_Database = -1;
                this.siteBase_Database = -1;

            }

            //Control
            foreach (SNP s in snps)
            {
                if (s.DatabaseParent != databaseParent || s.SiteBasis_Database != siteBase_Database || s.SiteIndelShift_Database != this.siteIndelShift_Database || s.IsFullSNP != this.snpAlleles
                    || this.databaseCharacter != s.DatabaseChar
                    ) throw new InvalidOperationException("All SNPs for a common loci must rever to the same database parent and database position - otherwise the loci would not be equal");
            }
            this.allAlleles = GetAlleleDictionary(all);
        }

        /// <summary>
        /// BASIC variables
        /// </summary>

        private int siteBase_Database;
        private int siteIndelShift_Database;
        private string databaseParent;

        /// <summary>
        /// Additional variables
        /// </summary>
        private List<SNP> all;
        private List<SNP> invalid;
        private Dictionary<char, List<SNP>> allAlleles;
        private bool snpAlleles;
        private char databaseCharacter;
        private int databaseParentSNPCount;

        private ISNPSiteFormater formater;


        private string GetSiteIDDatabase()
        {
            return String.Format("{0}-{1}", this.siteBase_Database, this.siteIndelShift_Database);
        }

        private void SetSiteIDDatabase(string snpSiteID_database)
        {
            string[] temp = snpSiteID_database.Split(new char[] { '-' });

            this.siteBase_Database = Convert.ToInt32(temp[0]);
            this.siteIndelShift_Database = Convert.ToInt32(temp[1]);
        }

        /// <summary>
        /// Get or set the position of this SNP site in the database sequence
        /// </summary>
        public string SiteID_Database
        {
            get
            {
                return this.GetSiteIDDatabase();
            }
            set
            {
                this.SetSiteIDDatabase(value);
            }
        }

        public int SiteBase_Database
        {
            get
            {
                return this.siteBase_Database;
            }
            set
            {
                this.siteBase_Database = value;
            }
        }

        public int SiteIndelShift_Database
        {
            get
            {
                return this.siteIndelShift_Database;
            }
            set
            {
                this.siteIndelShift_Database = value;
            }
        }



        /// <summary>
        /// Get a collection containing all valid SNPs
        /// </summary>
        public List<SNP> SNPs_Valid
        {
            get
            {
                return all;
            }
        }

        /// <summary>
        /// Get a collection containing all invalid SNPs
        /// </summary>
        public List<SNP> SNPs_Invalid
        {
            get
            {
                return this.invalid;
            }
        }


        /// <summary>
        /// Get or set the name of the database parent
        /// </summary>
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

        /// <summary>
        /// The number of SNPs which share this database parent.
        /// May be used to calculate the prominence of a locus.
        /// </summary>
        public int Number_SNPsAtDatabaseParent
        {
            get
            {
                return this.databaseParentSNPCount;
            }
            set
            {
                this.databaseParentSNPCount = value;
            }

        }

        /// <summary>
        /// Get the number of SNPs at this SNP site; only valid SNPs are considered
        /// </summary>
        public int Count_Valid
        {
            get
            {
                return this.all.Count;
            }

        }



        /// <summary>
        /// Indicates whether the SNPs at this SNP site represent fullSNPs or halfSNPs
        /// </summary>
        public bool IsFullSNP
        {
            get
            {
                return this.snpAlleles;
            }
        }




        /// <summary>
        /// Get the total number of SNP alleles found at this SNP site.
        /// </summary>
        public int Count_Alleles
        {
            get
            {

                return this.Alleles.Count;
            }
        }



        /// <summary>
        /// Get a dictionary containing all SNP-alleles at the given locus.
        /// The dictionary is calculated for all SNPs at this locus (contrast: Alleles_Subset).
        /// </summary>
        public Dictionary<char, List<SNP>> Alleles
        {
            get
            {
                return this.allAlleles;

            }

        }

        /// <summary>
        /// Get the number of all SNPs at this SNP site.
        /// </summary>
        public int Count_ValidAndInvalid
        {
            get
            {
                return this.all.Count + this.invalid.Count;
            }
        }

        /// <summary>
        /// Get the number of all invalid SNPs at this SNP site
        /// </summary>
        public int Count_Invalid
        {
            get
            {
                return this.invalid.Count;
            }
        }

        /// <summary>
        /// Get the percentage of the SNPs being valid. 
        /// </summary>
        public double PercentValidSNPs
        {
            get
            {
                return 100.0 * (double)this.all.Count / (double)this.Count_ValidAndInvalid;
            }
        }

        /// <summary>
        /// Obtain an allele dictionary
        /// </summary>
        /// <param name="snps"></param>
        /// <returns></returns>
        private Dictionary<char, List<SNP>> GetAlleleDictionary(List<SNP> snps)
        {
            Dictionary<char, List<SNP>> dict = new Dictionary<char, List<SNP>>();
            if (snps.Count == 0) return dict;
            List<SNP> temp = new List<SNP>();
            snps.Sort(new Sort.SortSNPs_DatabaseName_DatabasePosition_Allele());

            temp.Add(snps[0]);
            for (int i = 1; i < snps.Count; i++)
            {
                SNP lastSnp = temp[temp.Count - 1];
                if (lastSnp.QueryChar == snps[i].QueryChar)
                {
                    temp.Add(snps[i]);
                }
                else
                {
                    dict.Add(lastSnp.QueryChar, temp);
                    temp = new List<SNP>();
                    temp.Add(snps[i]);
                }
            }
            dict.Add(temp[0].QueryChar, temp); //also add the last
            return dict;
        }


        /// <summary>
        /// Relative importance of this particular SNP site with respect to the number of all SNPs mapping to this database parent.
        /// </summary>
        public double SNPSiteProminence_Percent
        {
            get
            {
                if (this.all.Count == 0) return 0.0;
                else return 100.0 * (double)this.Count_ValidAndInvalid / (double)this.databaseParentSNPCount;
            }
        }

        /// <summary>
        /// Set the active SNP subset at the given locus
        /// </summary>
        /// <param name="showPlusPlus">consider SNPs from the sense alignment</param>
        /// <param name="showPlusMinus">consider SNPs from the antisense alignment</param>
        /// <param name="show_A">consider SNPs having the nucleotide 'A'</param>
        /// <param name="show_T">consider SNPs having the nucleotide 'T'</param>
        /// <param name="show_C">consider SNPs having the nucleotide 'C'</param>
        /// <param name="show_G">consider SNPs having the nucleotide 'G'</param>
        public SNPSite GetSubset(bool showPlusPlus, bool showPlusMinus, bool show_A, bool show_T, bool show_C, bool show_G, bool show_indel, int minimumDistanceFromAlignmentEnd, int? maxLowQualityTokens, int? minimumSequenceQualitySite, int? minimumSequenceQualityNeighborhood)
        {
            List<SNP> subset = new List<SNP>();
            foreach (SNP s in this.all)
            {
                if (            //(s.DistanceFromAlignmentEnd==null) || no SNP should have a distance from the alignment end of null!!
                    (s.DistanceFromAlignmentEnd.Value >= minimumDistanceFromAlignmentEnd)) //Check distance
                {
                    if ((showPlusPlus && s.PlusPlusStrand) || (showPlusMinus && !s.PlusPlusStrand)) //Check strand
                    {

                        if ((maxLowQualityTokens == null) || (s.AlignmentQualityTokens <= maxLowQualityTokens.Value))
                        {
                            if ((minimumSequenceQualitySite == null || (s.SequenceQualityAtSite != null && s.SequenceQualityAtSite >= minimumSequenceQualitySite))
                                && (minimumSequenceQualityNeighborhood == null || (s.SequenceQualityInNeighborhood != null && s.SequenceQualityInNeighborhood >= minimumSequenceQualityNeighborhood)))
                            {
                                switch (s.QueryChar) //Check character state
                                {
                                    case 'A':
                                    case 'a': if (show_A) subset.Add(s);
                                        break;
                                    case 'T':
                                    case 't': if (show_T) subset.Add(s);
                                        break;
                                    case 'C':
                                    case 'c': if (show_C) subset.Add(s);
                                        break;
                                    case 'G':
                                    case 'g': if (show_G) subset.Add(s);
                                        break;
                                    case '-': if (show_indel) subset.Add(s);
                                        break;
                                    default: throw new InvalidOperationException("This is not possible");
                                }
                            }
                        }

                    }
                }
            }
            return new SNPSite(subset);
        }

        /// <summary>
        /// Set the active SNP subset at the given locus
        /// </summary>
        /// <param name="showPlusPlus">consider SNPs from the sense alignment</param>
        /// <param name="showPlusMinus">consider SNPs from the antisense alignment</param>
        /// <param name="show_A">consider SNPs having the nucleotide 'A'</param>
        /// <param name="show_T">consider SNPs having the nucleotide 'T'</param>
        /// <param name="show_C">consider SNPs having the nucleotide 'C'</param>
        /// <param name="show_G">consider SNPs having the nucleotide 'G'</param>
        public SNPSite GetSubset(bool showPlusPlus, bool showPlusMinus, bool show_A, bool show_T, bool show_C, bool show_G, bool show_indel, int minimuDistanceFromAlignmentEnd, int? maxLowQualityTokens, int? minimumSequenceQualitySite, int? minimumSequenceQualityNeighborhood, int minimumFrequencyOfAlleles, int showOnlyTheNBestAlleles)
        {
            SNPSite tempSite = GetSubset(showPlusPlus, showPlusMinus, show_A, show_T, show_C, show_G, show_indel, minimuDistanceFromAlignmentEnd, maxLowQualityTokens, minimumSequenceQualitySite, minimumSequenceQualityNeighborhood); //Set the subset with the already defined method
            Dictionary<char, List<SNP>> dict = tempSite.Alleles; //use the subset of SNP alleles 

            ////Create a subAllele dictionary containing only those alleles having at least the frequency minimumFrequncyOfAlleles
            List<List<SNP>> temp = new List<List<SNP>>();
            foreach (KeyValuePair<char, List<SNP>> kv in dict)
            {
                if (kv.Value.Count >= minimumFrequencyOfAlleles) temp.Add(kv.Value); //Minimum frequency of alleles is handled;
            }

            //Set the new suballeles according to the allele count
            temp.Sort(new Sort.SortSNPList_Count());

            //Create both a new subset of SNPs and a new subset of AlleleDictionary
            List<SNP> subset = new List<SNP>();
            int count = 0;
            foreach (List<SNP> list in temp)
            {
                if (count < showOnlyTheNBestAlleles)
                {
                    subset.AddRange(list);
                    count++;
                }
            }

            return new SNPSite(subset);
        }
        /// <summary>
        /// Set the active SNP subset at the given locus
        /// </summary>
        /// <param name="showPlusPlus">consider SNPs from the sense alignment</param>
        /// <param name="showPlusMinus">consider SNPs from the antisense alignment</param>
        /// <param name="show_A">consider SNPs having the nucleotide 'A'</param>
        /// <param name="show_T">consider SNPs having the nucleotide 'T'</param>
        /// <param name="show_C">consider SNPs having the nucleotide 'C'</param>
        /// <param name="show_G">consider SNPs having the nucleotide 'G'</param>
        public SNPSite GetSubset(bool showPlusPlus, bool showPlusMinus, bool show_A, bool show_T, bool show_C, bool show_G, bool showIndel, int minimuDistanceFromAlignmentEnd, int? maxLowQualityTokens, int? minimumSequenceQualitySite, int? minimumSequenceQualityNeighborhood, int minimumFrequencyOfAlleles, double minimumAllelFrequencyPercent, int showOnlyTheNBestAlleles)
        {
            double minFreqConverted = minimumAllelFrequencyPercent / 100.0;
            SNPSite tempSite = GetSubset(showPlusPlus, showPlusMinus, show_A, show_T, show_C, show_G, showIndel, minimuDistanceFromAlignmentEnd, maxLowQualityTokens, minimumSequenceQualitySite, minimumSequenceQualityNeighborhood); //Set the subset with the already defined method
            Dictionary<char, List<SNP>> dict = tempSite.Alleles; //use the subset of SNP alleles 

            ////Create a subAllele dictionary containing only those alleles having at least the frequency minimumFrequncyOfAlleles
            List<List<SNP>> temp = new List<List<SNP>>();
            foreach (KeyValuePair<char, List<SNP>> kv in dict)
            {
                double alleleFrequency = (double)kv.Value.Count / (double)this.Count_ValidAndInvalid;
                if (kv.Value.Count >= minimumFrequencyOfAlleles && alleleFrequency >= minFreqConverted) temp.Add(kv.Value); //Minimum frequency of alleles is handled;
            }

            //Set the new suballeles according to the allele count
            temp.Sort(new Sort.SortSNPList_Count());

            //Create both a new subset of SNPs and a new subset of AlleleDictionary
            List<SNP> subset = new List<SNP>();
            int count = 0;
            foreach (List<SNP> list in temp)
            {
                if (count < showOnlyTheNBestAlleles)
                {
                    subset.AddRange(list);
                    count++;
                }
            }

            return new SNPSite(subset);
        }


        /// <summary>
        /// Return the Polymorphism Information Content of all SNPs at this SNP site.
        /// </summary>
        public double PIC
        {
            get
            {


                double pic = 1.0;
                foreach (List<SNP> ls in this.Alleles.Values)
                {
                    pic -= Math.Pow((double)ls.Count / (double)this.Count_Valid, 2.0);
                }
                return pic;
            }
        }


        /// <summary>
        /// Get the number of transitions at this SNP site
        /// </summary>
        public int Count_Transitions
        {
            get
            {
                int count = 0;
                foreach (SNP s in all)
                {
                    if (s.Transition) count++;
                }
                return count;
            }
        }

        /// <summary>
        /// Get the number of transversions at this SNP site
        /// </summary>
        public int Count_Transversions
        {
            get
            {

                int count = 0;
                foreach (SNP s in all)
                {
                    if (s.Transversion) count++;
                }
                return count;

            }
        }

        /// <summary>
        /// Get the number of plus/plus SNPs at this SNP site
        /// </summary>
        public int Count_PlusPlus
        {
            get
            {
                int count = 0;
                foreach (SNP s in all)
                {
                    if (s.PlusPlusStrand) count++;
                }
                return count;
            }
        }


        /// <summary>
        /// Get the character state of the database sequence
        /// </summary>
        public char DatabaseCharacter
        {
            get
            {
                return this.databaseCharacter;
            }
        }

        /// <summary>
        /// Get the number of plus/minus SNPs at this SNP site
        /// </summary>
        public int Count_PlusMinus
        {
            get
            {
                int count = 0;
                foreach (SNP s in all)
                {
                    if (!s.PlusPlusStrand) count++;
                }
                return count;
            }
        }

        /// <summary>
        /// Get or set the StartSite of the SNP
        /// </summary>
        public int? Start
        {
            get
            {
                return this.SiteBase_Database;
            }
            set
            {
                this.SiteBase_Database = value.Value;
            }
        }

        public int? End
        {
            get
            {
                return this.SiteBase_Database;
            }
            set
            {
                this.SiteBase_Database = value.Value;
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
                if (value.Value != false) throw new InvalidOperationException("A SNP can never be the root!");
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

        /// <summary>
        /// Length of a SNP site is always 1
        /// </summary>
        public long Length
        {
            get { return 1; }

        }
        public ISNPSiteFormater Formater
        {
            get
            {

                return this.formater;
            }
            set
            {
                this.formater = value;
            }
        }


        public override string ToString()
        {
            if (this.formater != null) return formater.GetFormated(this);
            else return null;
        }





        public string FeatureName
        {
            get {return "SNP Site"; }
        }
    }



    public interface ISNPExtractor
    {
        List<SNP> GetSNPs(IPairwiseNucleotideSequenceAlignmentContainer pwc);
    }



    /// <summary>
    /// Extract SNPs from common nucleotide sequences
    /// </summary>
    public class HalfSNPExtractor : ISNPExtractor
    {




        public HalfSNPExtractor()
        {
        }



        /// <summary>
        /// Retrieve the SNPs for the given PairwiseNucleotideSequenceAlignment
        /// </summary>
        /// <param name="pw"></param>
        /// <returns></returns>
        private List<SNP> GetSNPsIntern(PairwiseNucleotideSequenceAlignment pw)
        {
            List<SNP> snps = new List<SNP>();
            //If the snps have not yet been extracted-> extract them

                ISequenceContainer data = pw.Alignment.DatabaseSequence;
                ISequenceContainer query = pw.Alignment.QuerySequence;

                DNASequenceValidator_ATCG val = new DNASequenceValidator_ATCG();
                int countData = -1;
                int countQuery = -1;
                List<SNP> toRet = new List<SNP>();

                for (int i = 0; i < data.Length; i++)
                {
                    if (data[i] != '-') countData++;
                    if (query[i] != '-') countQuery++;

                    //Test fundamental SNP characteristics
                    if (data[i] != query[i] && val.IsValid(data[i]) && val.IsValid(query[i]))
                    {
                        int minimumEndDistance = Math.Min(i, data.Length - i - 1);
                        SNP s = new SNP(pw.DatabaseParent, pw.QueryParent, data[i], query[i], pw.StartDatabase + countData, 0, pw.StartQuery + countQuery, 0, pw.PlusPlusStrand, false, true, minimumEndDistance);
                        toRet.Add(s);
                    }

                }
                snps = toRet;
            
            return snps;

        }





        public List<SNP> GetSNPs(IPairwiseNucleotideSequenceAlignmentContainer pwc)
        {
            List<SNP> toRet = new List<SNP>();
            List<IPairwiseAlignmentContainer> alignments = pwc.SubAlignments;
            foreach (IPairwiseAlignmentContainer pwac in alignments)
            {
                List<SNP> toAdd = new List<SNP>();

                //if the Pairwise alignment container is a leeve
                if (pwac is PairwiseNucleotideSequenceAlignment)
                {
                    toAdd = GetSNPsIntern((PairwiseNucleotideSequenceAlignment)pwac);
                }//if the container is another branch recursively call the method
                else if (pwac is IPairwiseNucleotideSequenceAlignmentContainer) toAdd = GetSNPs((IPairwiseNucleotideSequenceAlignmentContainer)pwac);
                else throw new InvalidOperationException("pwc has to be a PairiwiseNucleotide sequence container");

                //if the number of identified SNPs is larger than one add the SNPs
                if (toAdd.Count > 0) toRet.AddRange(toAdd);
            }

            return toRet;

        }

    }

    /// <summary>
    /// Extract SNPs from common nucleotide sequences
    /// </summary>
    public class HalfSNPExtractor_Indel : ISNPExtractor
    {




        public HalfSNPExtractor_Indel()
        {
        }



        /// <summary>
        /// Retrieve the SNPs for the given PairwiseNucleotideSequenceAlignment
        /// </summary>
        /// <param name="pw"></param>
        /// <returns></returns>
        private List<SNP> GetSNPsIntern(PairwiseNucleotideSequenceAlignment pw)
        {
            List<SNP> snps = new List<SNP>();
            //If the snps have not yet been extracted-> extract them

                ISequenceContainer data = pw.Alignment.DatabaseSequence;
                ISequenceContainer query = pw.Alignment.QuerySequence;

                DNASequenceValidator_ATCGIndel val = new DNASequenceValidator_ATCGIndel();
                int countData = -1;
                int countQuery = -1;
                int countIndelData = 0;
                int countIndelQuery = 0;
                List<SNP> toRet = new List<SNP>();

                for (int i = 0; i < data.Length; i++)
                {
                    if (data[i] != '-')
                    {
                        countData++;
                        countIndelData = 0;
                    }
                    else
                    {
                        countIndelData++;
                    }

                    if (query[i] != '-')
                    {
                        countQuery++;
                        countIndelQuery = 0;
                    }
                    else
                    {
                        countIndelQuery++;
                    }

                    //Test fundamental SNP characteristics
                    if (data[i] != query[i] && val.IsValid(data[i]) && val.IsValid(query[i]))
                    {
                        int minimumEndDistance = Math.Min(i, data.Length - i - 1);
                        SNP s = new SNP(pw.DatabaseParent, pw.QueryParent, data[i], query[i], pw.StartDatabase + countData, countIndelData, pw.StartQuery + countQuery, countIndelQuery, pw.PlusPlusStrand, false, true, minimumEndDistance);
                        toRet.Add(s);
                    }

                }
                snps = toRet;
            
            return snps;

        }


        public List<SNP> GetSNPs(IPairwiseNucleotideSequenceAlignmentContainer pwc)
        {
            List<SNP> toRet = new List<SNP>();
            List<IPairwiseAlignmentContainer> alignments = pwc.SubAlignments;
            foreach (IPairwiseAlignmentContainer pwac in alignments)
            {
                List<SNP> toAdd = new List<SNP>();

                //if the Pairwise alignment container is a leeve
                if (pwac is PairwiseNucleotideSequenceAlignment)
                {
                    toAdd = GetSNPsIntern((PairwiseNucleotideSequenceAlignment)pwac);
                }//if the container is another branch recursively call the method
                else if (pwac is IPairwiseNucleotideSequenceAlignmentContainer) toAdd = GetSNPs((IPairwiseNucleotideSequenceAlignmentContainer)pwac);
                else throw new InvalidOperationException("pwc has to be a PairiwiseNucleotide sequence container");

                //if the number of identified SNPs is larger than one add the SNPs
                if (toAdd.Count > 0) toRet.AddRange(toAdd);
            }

            return toRet;

        }

    }

    public class FullSNPExtractor
    {
        private int minDistanceFromEnd = 10;
        private int minAlleleFrequencyCount;
        private float minAlleleFrequencyPercent;
        private int minAlleles;
        private int maxAlleles;
        private int minCoverage;


        private List<IPairwiseNucleotideSequenceAlignmentContainer> pwa;
        private bool considerIndel;
        private ISNPQualityAlignment quality;
        private int countHalfSNPs = 0;
        private int countHalfSNPsites = 0;
        private int countAlignmentLength = 0;
        private int countInitialFullSNPs = 0;
        private int countInitialFullSNPSites = 0;
        private int countEndFullSNPs = 0;
        private int countEndFullSNPSites = 0;
        private List<SNP> toRet;


        public FullSNPExtractor(List<IPairwiseNucleotideSequenceAlignmentContainer> pwa, bool considerIndel, ISNPQualityAlignment quality, int minAllelFrequencyCount, float minAlleleFrequencyPercent, int minAlleles, int maxAlleles, int minCoverage)
        {
            this.considerIndel = considerIndel;
            this.pwa = pwa;
            this.quality = quality;
            this.minAlleleFrequencyCount = minAllelFrequencyCount;
            this.minAlleleFrequencyPercent = minAlleleFrequencyPercent;
            this.minAlleles = minAlleles;
            this.maxAlleles = maxAlleles;
            this.minCoverage = minCoverage;
            ExtractSNPs();

        }

        private void ExtractSNPs()
        {
            //StepOne HalfSNPs

            ISNPExtractor halfSNPExtractor;
            if (considerIndel) halfSNPExtractor = new HalfSNPExtractor_Indel();
            else halfSNPExtractor = new HalfSNPExtractor();
            List<SNP> halfSNPs = new List<SNP>();
            foreach (IPairwiseNucleotideSequenceAlignmentContainer pw in pwa)
            {
                countAlignmentLength += (int)pw.Length;
                halfSNPs.AddRange(halfSNPExtractor.GetSNPs(pw));
            }
            this.countHalfSNPs = halfSNPs.Count;


            //STEP 2 delimit  halfSNPs;
            List<SNPSite> halfSNPSites = SNPUtility.GetSNPSites(halfSNPs);
            this.countHalfSNPsites = halfSNPSites.Count;
            List<SNP> snpsForExtraction = new List<SNP>();
            foreach (SNPSite ss in halfSNPSites)
            {
                SNPSite tempSite = ss.GetSubset(true, true, true, true, true, true, true, minDistanceFromEnd, null, null, null);
                Dictionary<char, List<SNP>> alleles = tempSite.Alleles;
                bool fullfilledRequirement = false;

                foreach (List<SNP> sn in alleles.Values)
                {
                    if (sn.Count >= minAlleleFrequencyCount) fullfilledRequirement = true;
                }
                if (fullfilledRequirement)
                {
                    snpsForExtraction.AddRange(tempSite.SNPs_Valid);
                }
            }

            //STEP 3
            Dictionary<string, List<SNPSite>> locusDict = SNPUtility.GetSNPSiteDictionary_DatabaseName(snpsForExtraction);
            foreach (KeyValuePair<string, List<SNPSite>> kvp in locusDict)
            {
                this.countInitialFullSNPSites += kvp.Value.Count;
            }
            ISNPExtractor fullExtractor;
            if (considerIndel) fullExtractor = new FullSNPExtractorSite_Indel(locusDict);
            else fullExtractor = new FullSNPExtractorSite(locusDict);
            List<SNP> fullSNPs = new List<SNP>();
            foreach (IPairwiseNucleotideSequenceAlignmentContainer pw in pwa)
            {

                List<SNP> tempSNPs = fullExtractor.GetSNPs(pw);
                foreach (SNP sn in tempSNPs)
                {
                    if (quality != null) sn.AlignmentQualityTokens = this.quality.GetLowQualityToken(sn, pw);
                }
                fullSNPs.AddRange(tempSNPs);
            }
            this.countInitialFullSNPs = fullSNPs.Count;


            //STEP 4 Delimit the SNP-Alleles
            List<SNPSite> fullSites = SNPUtility.GetSNPSites(fullSNPs);
            toRet = new List<SNP>();

            foreach (SNPSite ss in fullSites)
            {
                SNPSite tempSite = ss.GetSubset(true, true, true, true, true, true, considerIndel, 0, null, null, null);

                Dictionary<char, List<SNP>> alleles = tempSite.Alleles;
                int allelesFullfillingRequirement = 0;

                foreach (List<SNP> sn in alleles.Values)
                {
                    float freqPercen = (float)100.0 * sn.Count / ss.Count_ValidAndInvalid;
                    if (sn.Count >= this.minAlleleFrequencyCount && freqPercen >= minAlleleFrequencyPercent) allelesFullfillingRequirement++;
                }

                if (allelesFullfillingRequirement >= minAlleles && allelesFullfillingRequirement <= maxAlleles && ss.Count_ValidAndInvalid >= minCoverage)
                {
                    toRet.AddRange(ss.SNPs_Valid);
                    toRet.AddRange(ss.SNPs_Invalid);
                    this.countEndFullSNPSites++;
                }
            }
            this.countEndFullSNPs = toRet.Count;
        }

        public List<SNP> GetSNPs
        {
            get
            {
                return this.toRet;
            }
        }

        public int Count_HalfSNPs
        {
            get
            {
                return countHalfSNPs;
            }
        }
        public int Count_HalfSNPSites
        {
            get
            {
                return countHalfSNPsites;
            }
        }
        public int Count_InitialFullSNPs
        {
            get
            {
                return countInitialFullSNPs;
            }
        }
        public int Count_InitialFullSNPSites
        {
            get
            {
                return countInitialFullSNPSites;
            }
        }
        public int Count_FinalFullSNPs
        {
            get
            {
                return countEndFullSNPs;
            }
        }
        public int Count_FinalFullSNPSites
        {
            get
            {
                return countEndFullSNPSites;
            }
        }

        public int Count_AlignmentLength
        {
            get
            {
                return countAlignmentLength;
            }
        }

    }

    /// <summary>
    /// Extract SNP-Alleles from a given collection of SNP sites; Indels will not be considered as valid alleles
    /// </summary>
    public class FullSNPExtractorSite : ISNPExtractor
    {
        private Dictionary<string, List<SNPSite>> snpSites;


        public FullSNPExtractorSite(Dictionary<string, List<SNPSite>> snpSites)
        {

            this.snpSites = snpSites;

        }


        private List<SNP> GetSNPsIntern(PairwiseNucleotideSequenceAlignment pw)
        {
            ///Test if every member fullfills the requirements;
            List<SNPSite> sites;
            List<SNP> ret = new List<SNP>();
            DNASequenceValidator_ATCG val = new DNASequenceValidator_ATCG();
            if (this.snpSites.ContainsKey(pw.DatabaseParent))
            {
                sites = this.snpSites[pw.DatabaseParent];
            }
            else return ret;


            ///Everything OK continue
            //Convert the list of snp sites into a dictionary for fast and easy retrieval of the appropriate values
            Dictionary<int, byte> dict = new Dictionary<int, byte>(sites.Count * 2);
            foreach (SNPSite s in sites)
            {

                //Add the SNP site only if the indel shift is zero, i.e no indels occure
                if (s.SiteIndelShift_Database == 0 && !dict.ContainsKey(s.SiteBase_Database))
                {
                    dict.Add(s.SiteBase_Database, (byte)'0');
                }
            }

            //Set working variables
            ISequenceContainer data = pw.Alignment.DatabaseSequence;
            ISequenceContainer query = pw.Alignment.QuerySequence;

            int databaseCount = -1;
            int queryCount = -1;
            int countQueryShift = 0;

            for (int i = 0; i < data.Length; i++)
            {

                if (query[i] != '-')
                {
                    queryCount++;
                    countQueryShift = 0;
                }
                else countQueryShift++;

                //If the dictionary contains the key extract the SNP
                if (data[i] != '-')
                {
                    databaseCount++;

                    if (dict.ContainsKey(pw.StartDatabase + databaseCount))
                    {

                        int minimumEndDistance = Math.Min(i, data.Length - 1 - i);
                        //Test whether or not the SNP at the given SNP-locus is valid, if not ad the SNP anyway and raise the invalidAllel flag;
                        if (val.IsValid(data[i]) && val.IsValid(query[i]))
                        {

#if DEBUG
                            if (countQueryShift != 0) throw new InvalidOperationException("The query shift should be zero");
#endif
                            SNP snp = new SNP(pw.DatabaseParent, pw.QueryParent, data[i], query[i], pw.StartDatabase + databaseCount, 0, pw.StartQuery + queryCount, countQueryShift, pw.PlusPlusStrand, true, true, minimumEndDistance);
                            ret.Add(snp);
                        }
                        else
                        {
                            //Test for query indel shift
                            SNP snp = new SNP(pw.DatabaseParent, pw.QueryParent, data[i], query[i], pw.StartDatabase + databaseCount, 0, pw.StartQuery + queryCount, countQueryShift, pw.PlusPlusStrand, true, false, minimumEndDistance);
                            ret.Add(snp);
                        }

                    }
                }

            }

            return ret;


        }


        public List<SNP> GetSNPs(IPairwiseNucleotideSequenceAlignmentContainer pwc)
        {
            List<SNP> toRet = new List<SNP>();
            List<IPairwiseAlignmentContainer> alignments = pwc.SubAlignments;
            foreach (IPairwiseAlignmentContainer pwca in alignments)
            {
                List<SNP> toAdd = new List<SNP>();

                //if the Pairwise alignment container is a leeve
                if (pwca is PairwiseNucleotideSequenceAlignment)
                {
                    toAdd = GetSNPsIntern((PairwiseNucleotideSequenceAlignment)pwca);
                }//if the container is another branch recursively call the method
                else if (pwca is IPairwiseNucleotideSequenceAlignmentContainer) toAdd = GetSNPs((IPairwiseNucleotideSequenceAlignmentContainer)pwca);
                else throw new InvalidOperationException("pwc has to be a PairiwiseNucleotide sequence container");

                //if the number of identified SNPs is larger than one add the SNPs
                if (toAdd.Count > 0) toRet.AddRange(toAdd);
            }

            return toRet;

        }


    }



    /// <summary>
    /// Extract SNP-Alleles from a given collection of SNP sites; Indels will be considered as valid alleles
    /// </summary>
    public class FullSNPExtractorSite_Indel : ISNPExtractor
    {
        private Dictionary<string, List<SNPSite>> snpSites;

        public FullSNPExtractorSite_Indel(Dictionary<string, List<SNPSite>> snpSites)
        {

            this.snpSites = snpSites;
        }

        private List<SNP> GetSNPsIntern(PairwiseNucleotideSequenceAlignment pw)
        {
            ///Test if every member fullfills the requirements;
            List<SNPSite> sites;
            List<SNP> ret = new List<SNP>();
            DNASequenceValidator_ATCGIndel val = new DNASequenceValidator_ATCGIndel();
            if (this.snpSites.ContainsKey(pw.DatabaseParent))
            {
                sites = this.snpSites[pw.DatabaseParent];
            }
            else return ret;


            ///Everything OK continue
            //Convert the list of snp sites into a dictionary for fast and easy retrieval of the appropriate values
            Dictionary<int, Dictionary<int, bool>> dict = new Dictionary<int, Dictionary<int, bool>>(sites.Count * 2);
            foreach (SNPSite s in sites)
            {
                if (dict.ContainsKey(s.SiteBase_Database))
                {
                    dict[s.SiteBase_Database].Add(s.SiteIndelShift_Database, false);
                }
                else
                {
                    Dictionary<int, bool> d = new Dictionary<int, bool>();
                    d.Add(s.SiteIndelShift_Database, false);
                    dict.Add(s.SiteBase_Database, d);
                }
            }

            //Set working variables
            ISequenceContainer data = pw.Alignment.DatabaseSequence;
            ISequenceContainer query = pw.Alignment.QuerySequence;

            int databaseCount = -1;
            int queryCount = -1;
            int countQueryShift = 0;
            int countDataShift = 0;


            Dictionary<int, bool> activeDictionary = new Dictionary<int, bool>();
            int i = 0;
            for (; i < data.Length; i++)
            {





                //If the dictionary contains the key extract the SNP
                if (data[i] != '-')
                {
                    foreach (KeyValuePair<int, bool> kv in activeDictionary)
                    {
                        if (kv.Value == false)
                        {
                            int minimumEndDistance = Math.Min(i, data.Length - 1 - i);
                            int addShift = kv.Key - countDataShift;
                            SNP snp = new SNP(pw.DatabaseParent, pw.QueryParent, '-', '-', pw.StartDatabase + databaseCount, kv.Key, pw.StartQuery + queryCount, addShift, pw.PlusPlusStrand, true, true, minimumEndDistance);
                            ret.Add(snp);
                        }
                    }
                    databaseCount++;
                    countDataShift = 0;
                    if (i < data.Length && dict.ContainsKey(pw.StartDatabase + databaseCount))
                    {
                        activeDictionary = dict[pw.StartDatabase + databaseCount];
                    }
                    else activeDictionary = new Dictionary<int, bool>();

                }
                else
                {
                    countDataShift++;
                }

                if (query[i] != '-')
                {
                    queryCount++;
                    countQueryShift = 0;
                }
                else countQueryShift++;

                if (activeDictionary.ContainsKey(countDataShift))
                {
                    activeDictionary[countDataShift] = true;
                    int minimumEndDistance = Math.Min(i, data.Length - 1 - i);


                    //Test whether or not the SNP at the given SNP-locus is valid, if not ad the SNP anyway and raise the invalidAllel flag;
                    if (val.IsValid(data[i]) && val.IsValid(query[i]))
                    {
                        SNP snp = new SNP(pw.DatabaseParent, pw.QueryParent, data[i], query[i], pw.StartDatabase + databaseCount, countDataShift, pw.StartQuery + queryCount, countQueryShift, pw.PlusPlusStrand, true, true, minimumEndDistance);
                        ret.Add(snp);
                    }
                    else
                    {
                        SNP snp = new SNP(pw.DatabaseParent, pw.QueryParent, data[i], query[i], pw.StartDatabase + databaseCount, countDataShift, pw.StartQuery + queryCount, countQueryShift, pw.PlusPlusStrand, true, false, minimumEndDistance);
                        ret.Add(snp);
                    }
                }


            }



            ///Also process the last ones
            foreach (KeyValuePair<int, bool> kv in activeDictionary)
            {
                if (kv.Value == false)
                {
                    int minimumEndDistance = Math.Min(i, data.Length - 1 - i);
                    if (minimumEndDistance < 0) minimumEndDistance = 0;
                    int addShift = kv.Key - countDataShift;
                    SNP snp = new SNP(pw.DatabaseParent, pw.QueryParent, '-', '-', pw.StartDatabase + databaseCount, kv.Key, pw.StartQuery + queryCount, addShift, pw.PlusPlusStrand, true, true, minimumEndDistance);
                    ret.Add(snp);
                }
            }


            return ret;


        }

        public List<SNP> GetSNPs(IPairwiseNucleotideSequenceAlignmentContainer pwc)
        {
            List<SNP> toRet = new List<SNP>();
            List<IPairwiseAlignmentContainer> alignments = pwc.SubAlignments;
            foreach (IPairwiseAlignmentContainer pwca in alignments)
            {
                List<SNP> toAdd = new List<SNP>();

                //if the Pairwise alignment container is a leeve
                if (pwca is PairwiseNucleotideSequenceAlignment)
                {
                    toAdd = GetSNPsIntern((PairwiseNucleotideSequenceAlignment)pwca);
                }//if the container is another branch recursively call the method
                else if (pwca is IPairwiseNucleotideSequenceAlignmentContainer) toAdd = GetSNPs((IPairwiseNucleotideSequenceAlignmentContainer)pwca);
                else throw new InvalidOperationException("pwc has to be a PairiwiseNucleotide sequence container");

                //if the number of identified SNPs is larger than one add the SNPs
                if (toAdd.Count > 0) toRet.AddRange(toAdd);
            }

            return toRet;

        }


    }

    public interface ISNPSiteFormater
    {
        string GetFormated(SNPSite ss);
        string GetHeader();
    }

    public interface ISibblingSNPFormater
    {
        string GetFormated(SibblingSNPs ss);
        string GetHeader();
    }

    /// <summary>
    /// Interface for SNPs
    /// </summary>
    public interface ISNP : IPositionable
    {
        string DatabaseParent { get;set;}
        char DatabaseChar { get;}
        char QueryChar { get;}
        string SiteID_Database { get;set;}
        string SiteID_Query { get;set;}
        int SiteIndelShift_Database { get;set;}
        int SiteIndelShift_Query { get;set;}
        int SiteBasis_Database { get;set;}
        int SiteBasis_Query { get;set;}
        bool Transition { get;}
        bool Transversion { get;}
        string QueryParent { get;set;}
        IPairwiseNucleotideSequenceAlignmentContainer ParentAlignment { get;set;}
        bool PlusPlusStrand { get;set;}

        //SiteID
        //SiteIndelShift_Database
        //

    }

    public interface ISNPQualityAlignment
    {
        int GetLowQualityToken(SNP snp, IPairwiseNucleotideSequenceAlignmentContainer pwa);
    }

    /// <summary>
    /// Estimates the quality of SNPs by neighborhood to low quality characters, like 'N' or '-' 
    /// </summary>
    public class SNPAlignmentQualityAssessment : ISNPQualityAlignment
    {

        private int slidingWindowSize = 0;
        private ISequenceValidator val = new DNASequenceValidator_ATCG();

        public SNPAlignmentQualityAssessment(int slidingWindowSize)
        {
            this.slidingWindowSize = slidingWindowSize;
        }

        public SNPAlignmentQualityAssessment(int slidingWindowSize, ISequenceValidator validator)
        {
            this.slidingWindowSize = slidingWindowSize;
            this.val = validator;
        }


        public int GetLowQualityToken(SNP snp, IPairwiseNucleotideSequenceAlignmentContainer ipwa)
        {
            //Check if sequences are identical
            if (snp.DatabaseParent != ipwa.DatabaseParent || snp.QueryParent != ipwa.QueryParent) throw new InvalidOperationException("Parent sequences of the SNP and of the alignment have to be identical");

            //Retrieve sequence
            PairwiseNucleotideSequenceAlignment pwa = ipwa.GetAlignmentCoveringDatabasePosition(snp.SiteBasis_Database);
            if (pwa == null) throw new Exception("Could not find SNP position");
            ISequenceContainer data = pwa.Alignment.DatabaseSequence;
            ISequenceContainer query = pwa.Alignment.QuerySequence;

            //Convert snp position into alignment frame
            int countDatabase = -1;
            int countDatabaseShift = 0;
            int snpPositionInAlignment = 0;

            bool foundBasis = false;
            for (int i = 0; i < data.Length; i++)
            {
                if (data[i] == '-')
                {
                    countDatabaseShift++;
                }
                else if (!foundBasis)
                {
                    countDatabase++;
                    countDatabaseShift = 0;
                }
                else
                {
                    snpPositionInAlignment = i - 1;
                    break;
                }

                //
                //       
                //0123456789
                if (pwa.StartDatabase + countDatabase == snp.SiteBasis_Database) foundBasis = true;

                if (foundBasis && snp.SiteIndelShift_Database == countDatabaseShift)
                {
                    snpPositionInAlignment = i;
                    break;
                }
            }

            int overheadStart = slidingWindowSize / 2;
            int overheadEnd = slidingWindowSize - overheadStart;


            //Determine start position of search
            //
            //0123456789
            //  --S--
            int startPos = snpPositionInAlignment - overheadStart;
            if (startPos < 0) startPos = 0;

            //Determine end position of search
            int endPos = snpPositionInAlignment + overheadEnd - 1;
            if (endPos > data.Length - 1) endPos = data.Length - 1;

            int countLowQualityCharacters = 0;

            //SNP at 6 - Windowsize 7 startPos= 3 endPos=9
            //
            //0123456789012345
            //      S
            for (int i = startPos; i <= endPos; i++)
            {
                if (!val.IsValid(data[i]) || !val.IsValid(query[i])) countLowQualityCharacters++;
            }

            return countLowQualityCharacters;
        }


        public int SlidingWindowSize
        {
            get
            {
                return this.slidingWindowSize;
            }
        }

    }


    public interface ISNPSequenceQualityAssessment
    {
        int GetQualityAtNeighborhood(Bio.SNP.SNP snp, Bio.Seq.QualitySequence qs);
        int GetQualityAtSite(Bio.SNP.SNP snp, Bio.Seq.QualitySequence qs);
    }

    /// <summary>
    /// Estimates the quality of SNPs by neighborhood to low quality characters, like 'N' or '-' 
    /// </summary>
    public class SNPSequenceQualityAssessment : ISNPSequenceQualityAssessment
    {
        private int snpNeighborhood;

        public SNPSequenceQualityAssessment(int snpNeighborhood)
        {
            if (snpNeighborhood < 3) throw new InvalidOperationException("SNP neighborhood must be at least 3 bp");
            this.snpNeighborhood = snpNeighborhood;
        }

        public SNPSequenceQualityAssessment()
        {
            this.snpNeighborhood = 5;
        }

        /// <summary>
        /// Get the quality at a SNP position;
        /// Automatically considers whether the SNP was found on the plus/plus or plus/minus orientation of the sequence; This is important because BLASTN performs a reverse complement and measures the distance always from the 3' end
        /// </summary>
        /// <param name="snp"></param>
        /// <param name="qs"></param>
        /// <returns></returns>
        public int GetQualityAtSite(SNP snp, QualitySequence qs)
        {
            if (snp.QueryParent != qs.Name) throw new InvalidOperationException("Parent sequences of the SNP and of the alignment have to be identical");

            int leng = (int)qs.Length;
            int snpPosition;
            int toRet = 0;

            if (snp.PlusPlusStrand)
            {
                snpPosition = snp.SiteBasis_Query - 1;
                toRet = qs.Sequence[snpPosition];
                if (snp.SiteIndelShift_Query > 0 && snpPosition + 1 < qs.Length)
                {
                    toRet += qs.Sequence[snpPosition + 1];
                    toRet = toRet / 2;
                }
            }
            else
            {   //     S
                //0123456789
                //1234567890 SNP with position positiv 6 is in reality at position 5 ie pos-1
                //0987654321 SNP with position negativ 2 is in reality at position 8 ie. length-pos
                snpPosition = leng - snp.SiteBasis_Query;
                toRet = qs.Sequence[snpPosition];
                if (snp.SiteIndelShift_Query > 0 && snpPosition - 1 >= 0)
                {
                    toRet += qs.Sequence[snpPosition - 1];
                    toRet = toRet / 2;
                }
            }

            return toRet;
        }



        public int GetQualityAtNeighborhood(SNP snp, QualitySequence qs)
        {
            //Check if sequences are identical
            if (snp.QueryParent != qs.Name) throw new InvalidOperationException("Parent sequences of the SNP and of the alignment have to be identical");
            ISequenceContainer qalS = qs.Sequence;

            int leng = qalS.Length;
            int snpPosition;
            if (snp.PlusPlusStrand)
            {
                snpPosition = snp.SiteBasis_Query - 1;
            }
            else
            {   //     S
                //0123456789
                //1234567890 SNP with position positiv 6 is in reality at position 5 ie pos-1
                //0987654321 SNP with position negativ 2 is in reality at position 8 ie. length-pos
                snpPosition = leng - snp.SiteBasis_Query;
            }

            int work = 0;
            int overheadStart = 0;
            int overheadEnd = 0;
            if (snp.SiteIndelShift_Query == 0)
            {
                work = snpNeighborhood;
                overheadStart = work / 2;
                overheadEnd = work - overheadStart;
            }
            else
            {
                work = snpNeighborhood - 1;
                if (snp.PlusPlusStrand)
                {

                    overheadStart = (work / 2) - 1;
                    overheadEnd = work - overheadStart;
                }
                else
                {
                    overheadStart = (work / 2);
                    overheadEnd = work - overheadStart;
                }
            }






            //Determine start position of search
            //
            //0123456789
            //  --S--
            int startPos = snpPosition - overheadStart;
            if (startPos < 0) startPos = 0;

            //Determine end position of search
            int endPos = snpPosition + overheadEnd - 1;
            if (endPos > qalS.Length - 1) endPos = qalS.Length - 1;

            int countSequenceQuality = 0;
            int countLeng = 0;

            //SNP at 6 - Windowsize 7 startPos= 3 endPos=9
            //
            //0123456789012345
            //      S
            for (int i = startPos; i <= endPos; i++)
            {

                if (i != snpPosition || snp.SiteIndelShift_Query > 0)
                {
                    countSequenceQuality += qalS[i];
                    countLeng++;
                }
            }

            return countSequenceQuality / countLeng;
        }

    }


    /// <summary>
    /// Estimates the quality of SNPs by neighborhood to low quality characters, like 'N' or '-' 
    /// </summary>
    public class SNPAlignmentQualityAssessment_454 : ISNPQualityAlignment
    {

        private int slidingWindowSize = 0;
        private int immidateNeighborhood = 3;
        private int penaltyNormal = 1;
        private int penaltyEqualCharacter = 3;
        private ISequenceValidator val = new DNASequenceValidator_ATCG();
        private SNPAlignmentQualityAssessment qualityAssessment;

        public SNPAlignmentQualityAssessment_454(int slidingWindowSize)
        {
            this.slidingWindowSize = slidingWindowSize;
            this.qualityAssessment = new SNPAlignmentQualityAssessment(slidingWindowSize);
        }

        public SNPAlignmentQualityAssessment_454(int slidingWindowSize, ISequenceValidator validator)
        {
            this.slidingWindowSize = slidingWindowSize;
            this.qualityAssessment = new SNPAlignmentQualityAssessment(slidingWindowSize, validator);
            this.val = validator;
        }


        public int GetLowQualityToken(SNP snp, IPairwiseNucleotideSequenceAlignmentContainer ipwa)
        {
            int countLowQualityCharacters = qualityAssessment.GetLowQualityToken(snp, ipwa);



            //Retrieve sequence
            PairwiseNucleotideSequenceAlignment pwa = ipwa.GetAlignmentCoveringDatabasePosition(snp.SiteBasis_Database);
            if (pwa == null) throw new Exception("Could not find SNP position");
            ISequenceContainer data = pwa.Alignment.DatabaseSequence;
            ISequenceContainer query = pwa.Alignment.QuerySequence;

            //Convert snp position into alignment frame
            int countDatabase = -1;
            int countShift = 0;
            int snpPositionInAlignment = 0;

            bool foundBasis = false;
            for (int k = 0; k < data.Length; k++)
            {
                if (data[k] == '-')
                {
                    countShift++;
                }
                else if (!foundBasis)
                {
                    countDatabase++;
                    countShift = 0;
                }
                else
                {
                    snpPositionInAlignment = k - 1;
                    break;
                }

                //
                //       
                //0123456789
                if (pwa.StartDatabase + countDatabase == snp.SiteBasis_Database) foundBasis = true;

                if (foundBasis && snp.SiteIndelShift_Database == countShift)
                {
                    snpPositionInAlignment = k;
                    break;
                }
            }

            int polyNpenalty = 0;
            char lastChar;




            // SNP at 2
            //0123456789 Length 10
            //  S
            int i = snpPositionInAlignment;

            while ((i - 1) >= 0 && snpPositionInAlignment - i <= this.immidateNeighborhood)
            {
                lastChar = query[i];

                if (lastChar != '-')
                {
                    while ((i - 1) >= 0 && (query[i - 1] == '-' || query[i - 1] == lastChar))
                    {
                        if (query[i - 1] != '-')
                        {
                            if (query[i - 1] == snp.QueryChar) polyNpenalty += penaltyEqualCharacter;
                            else polyNpenalty += penaltyNormal;
                        }
                        i--;
                    }
                }
                i--;
            }

            i = snpPositionInAlignment;

            while ((i + 1) <= query.Length - 1 && i - snpPositionInAlignment <= this.immidateNeighborhood)
            {
                lastChar = query[i];
                if (lastChar != '-')
                {
                    while ((i + 1) <= query.Length - 1 && (query[i + 1] == '-' || query[i + 1] == lastChar))
                    {
                        if (query[i + 1] != '-')
                        {
                            if (query[i + 1] == snp.QueryChar) polyNpenalty += penaltyEqualCharacter;
                            else polyNpenalty += penaltyNormal;
                        }
                        i++;
                    }
                }
                i++;
            }


            //Return the value
            return countLowQualityCharacters + polyNpenalty;
        }

        public int SlidingWindowSize
        {
            get
            {
                return this.slidingWindowSize;
            }
        }

    }


    /// <summary>
    /// A single nucleotide polymorphism
    /// Can either harbor a true SNP eg C<->A or an allelic SNP eg. A<->A
    /// </summary>
    public class SNP : ISNP
    {
        //BASIC characters
        private char databaseChar;
        private char queryChar;
        private string databaseName;
        private string queryName;
        private bool isSNPAllele = false;
        private bool valid = false;

        //Basic characters position
        private int siteIndelShift_database;
        private int siteIndelShift_query;
        private int siteBasis_database;
        private int siteBasis_query;


        //Additionals
        private object tag;
        private IPairwiseNucleotideSequenceAlignmentContainer parentAlignment;
        private bool? plusplusStrand = null;
        private string alleleID;

        //Quality related
        private int? distanceFromAlignmentEnd = null;
        private int? alignmentQualityToken = null;
        private int? sequenceQuality_queryPosition = null;
        private int? sequenceQuality_neighborhood = null;


        public SNP(string databaseParent, string queryParent, char databaseCharacter, char queryCharacter, string siteID_database, string siteID_query)
        {
            this.databaseName = databaseParent;
            this.queryName = queryParent;
            this.queryChar = queryCharacter;
            this.databaseChar = databaseCharacter;

            this.SetSiteIDDatabase(siteID_database);
            this.SetSiteIDQuery(siteID_query);
        }

        public SNP(string databaseParent, string queryParent, char databaseCharacter, char queryCharacter, string siteID_database, string siteID_query, bool plusPlusStrand)
            : this(databaseParent, queryParent, databaseCharacter, queryCharacter, siteID_database, siteID_query)
        {
            this.plusplusStrand = plusPlusStrand;
        }
        public SNP(string databaseParent, string queryParent, char databaseCharacter, char queryCharacter, string siteID_database, string siteID_query, bool plusPlusStrand, bool isFullSNP, bool isValid)
            : this(databaseParent, queryParent, databaseCharacter, queryCharacter, siteID_database, siteID_query, plusPlusStrand)
        {
            this.isSNPAllele = isFullSNP;
            this.valid = isValid;
        }

        public SNP(string databaseParent, string queryParent, char databaseCharacter, char queryCharacter, string siteID_database, string siteID_query, bool plusPlusStrand, bool isFullSNP, bool isValid, int distanceFromAlignmentEnd)
            : this(databaseParent, queryParent, databaseCharacter, queryCharacter, siteID_database, siteID_query, plusPlusStrand)
        {
            this.isSNPAllele = isFullSNP;
            this.valid = isValid;
            this.distanceFromAlignmentEnd = distanceFromAlignmentEnd;
        }

        public SNP(string databaseParent, string queryParent, char databaseCharacter, char queryCharacter, int siteBasis_database, int siteIndelShift_database, int siteBasis_query, int siteIndelShift_Query, bool plusPlusStrand, bool isFullSNP, bool isValid, int distanceFromAlignmentEnd)
        {
            this.databaseName = databaseParent;
            this.queryName = queryParent;
            this.databaseChar = databaseCharacter;
            this.queryChar = queryCharacter;
            this.siteBasis_database = siteBasis_database;
            this.siteIndelShift_database = siteIndelShift_database;
            this.siteBasis_query = siteBasis_query;
            this.siteIndelShift_query = siteIndelShift_Query;
            this.isSNPAllele = isFullSNP;
            this.valid = isValid;
            this.plusplusStrand = plusPlusStrand;
            this.distanceFromAlignmentEnd = distanceFromAlignmentEnd;
        }

        public SNP(string databaseParent, string queryParent, char databaseCharacter, char queryCharacter, string siteID_database, string siteID_query, bool plusPlusStrand, bool isFullSNP, bool isValid, int distanceFromAlignmentEnd, int lowAlignmentQualityToken, int sequenceQualitySite, int sequenceQualityNeighborhood)
            : this(databaseParent, queryParent, databaseCharacter, queryCharacter, siteID_database, siteID_query, plusPlusStrand)
        {
            this.isSNPAllele = isFullSNP;
            this.valid = isValid;
            this.distanceFromAlignmentEnd = distanceFromAlignmentEnd;
            this.alignmentQualityToken = lowAlignmentQualityToken;
            this.sequenceQuality_queryPosition = sequenceQualitySite;
            this.sequenceQuality_neighborhood = sequenceQualityNeighborhood;
        }

        public SNP(string databaseParent, string queryParent, char databaseCharacter, char queryCharacter, int siteBasis_database, int siteIndelShift_database, int siteBasis_query, int siteIndelShift_Query, bool plusPlusStrand, bool isFullSNP, bool isValid, int distanceFromAlignmentEnd, int lowAlignmentQualityToken, int sequenceQualitySite, int sequenceQualityNeighborhood)
            : this(databaseParent, queryParent, databaseCharacter, queryCharacter, siteBasis_database, siteIndelShift_database, siteBasis_query, siteIndelShift_Query, plusPlusStrand, isFullSNP, isValid, distanceFromAlignmentEnd)
        {
            this.alignmentQualityToken = lowAlignmentQualityToken;
            this.sequenceQuality_queryPosition = sequenceQualitySite;
            this.sequenceQuality_neighborhood = sequenceQualityNeighborhood;
        }

        public int? Start
        {
            get
            {
                return this.SiteBasis_Database;
            }
            set
            {
                this.SiteBasis_Database = value.Value;
            }
        }




        public int? End
        {
            get
            {
                return this.SiteBasis_Database;
            }
            set
            {
                throw new InvalidOperationException("It is not possible to set the end position of a SNP");
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
                if (value == true) throw new InvalidOperationException("A SNP cannot possible be the root, this value will therefore always be false");
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

        public int? SequenceQualityAtSite
        {
            get
            {
                return this.sequenceQuality_queryPosition;
            }
            set
            {
                this.sequenceQuality_queryPosition = value;
            }
        }

        public int? SequenceQualityInNeighborhood
        {
            get
            {
                return this.sequenceQuality_neighborhood;
            }
            set
            {
                this.sequenceQuality_neighborhood = value;
            }
        }

        public long Length
        {
            get { return 1; }
        }




        public string DatabaseParent
        {
            get
            {
                return this.databaseName;
            }
            set
            {
                this.databaseName = value;
            }
        }

        public char DatabaseChar
        {
            get
            {
                return this.databaseChar;
            }
        }

        public char QueryChar
        {
            get
            {
                return this.queryChar;
            }
        }

        /// <summary>
        /// the ID of the SNP site in the form 1221-2
        /// the first part is the basis and the second the indel shift
        /// </summary>
        public string SiteID_Database
        {
            get
            {

                return this.GetSiteIDDatabase();
            }
            set
            {
                this.SetSiteIDDatabase(value);
            }
        }

        public string SiteID_Query
        {
            get
            {

                return this.GetSiteIDQuery();
            }
            set
            {

                this.SetSiteIDQuery(value);
            }
        }

        /// <summary>
        /// Allele identification
        /// </summary>
        public string AlleleID
        {
            get
            {
                return this.alleleID;
            }
            set
            {
                this.alleleID = value;
            }
        }



        /// <summary>
        /// Indicates whether the SNP represents a transition.
        /// (Allelic SNPs are considered)
        /// </summary>
        public bool Transition
        {
            get
            {
                switch (this.databaseChar)
                {
                    case 'A':
                    case 'a': if (queryChar == 'g' || queryChar == 'G') return true;
                        break;
                    case 'T':
                    case 't': if (queryChar == 'C' || queryChar == 'c') return true;
                        break;
                    case 'G':
                    case 'g': if (queryChar == 'a' || queryChar == 'A') return true;
                        break;
                    case 'c':
                    case 'C': if (queryChar == 't' || queryChar == 'T') return true;
                        break;
                }
                return false;
            }

        }

        private string GetSiteIDDatabase()
        {
            return String.Format("{0}-{1}", this.siteBasis_database, this.siteIndelShift_database);
        }

        private void SetSiteIDDatabase(string snpSiteID_database)
        {
            string[] temp = snpSiteID_database.Split(new char[] { '-' });

            this.siteBasis_database = Convert.ToInt32(temp[0]);
            this.siteIndelShift_database = Convert.ToInt32(temp[1]);
        }

        private string GetSiteIDQuery()
        {
            return String.Format("{0}-{1}", this.siteBasis_query, this.siteIndelShift_query);
        }



        private void SetSiteIDQuery(string snpSiteID_query)
        {
            string[] temp = snpSiteID_query.Split(new char[] { '-' });

            this.siteBasis_query = Convert.ToInt32(temp[0]);
            this.siteIndelShift_query = Convert.ToInt32(temp[1]);
        }

        /// <summary>
        /// Indicates whether the SNP represents a transversion.
        /// (Allelic SNPs are considered)
        /// </summary>
        public bool Transversion
        {
            get
            {
                if (databaseChar == queryChar) return false; //For allelic snps
                return !this.Transition;
            }
        }

        public string QueryParent
        {
            get
            {
                return this.queryName;
            }
            set
            {
                this.queryName = value;
            }
        }

        public IPairwiseNucleotideSequenceAlignmentContainer ParentAlignment
        {
            get
            {
                return this.parentAlignment;
            }
            set
            {
                this.parentAlignment = value;
            }
        }

        /// <summary>
        /// Indicates whether the two strands are of the plus/plus alignment
        /// </summary>
        public bool PlusPlusStrand
        {
            get
            {
                return this.plusplusStrand.Value;
            }
            set
            {
                this.plusplusStrand = value;
            }
        }



        /// <summary>
        /// Indicates whether an SNP has been identified in course of an halfSNP or a fullSNP method;
        /// The fullSNPs are the biological meaningful - halfSNPs are just polymorphism between the database sequence and the query sequence
        /// </summary>
        public bool IsFullSNP
        {
            get
            {
                return this.isSNPAllele;
            }
            set
            {
                this.isSNPAllele = value;
            }
        }

        /// <summary>
        /// Indicates whether this SNP represents a valid allele; This should only be possible when the SNP-Allele identification algorithm was used
        /// </summary>
        public bool Valid
        {
            get
            {
                return this.valid;
            }
            set
            {
                this.valid = value;
            }
        }

        public int? DistanceFromAlignmentEnd
        {
            get
            {
                return this.distanceFromAlignmentEnd;
            }
            set
            {
                this.distanceFromAlignmentEnd = value;
            }
        }

        /// <summary>
        /// Get or set the estimated quality of the SNP, the higher the value the lower the quality
        /// </summary>
        public int? AlignmentQualityTokens
        {
            get
            {
                return this.alignmentQualityToken;
            }
            set
            {
                this.alignmentQualityToken = value;
            }
        }






        public int SiteIndelShift_Database
        {
            get
            {

                return this.siteIndelShift_database;
            }
            set
            {
                this.siteIndelShift_database = value;
            }
        }

        public int SiteIndelShift_Query
        {
            get
            {

                return this.siteIndelShift_query;
            }
            set
            {
                this.siteIndelShift_query = value;
            }
        }

        public int SiteBasis_Database
        {
            get
            {

                return this.siteBasis_database;
            }
            set
            {

                this.siteBasis_database = value;
            }
        }

        public int SiteBasis_Query
        {
            get
            {
                return this.siteBasis_query;
            }
            set
            {
                this.siteBasis_query = value;
            }
        }



        #region IFeature Members


        public string FeatureName
        {
            get { return "SNP"; }
        }

        #endregion
    }


    /// <summary>
    ///Contains a collection of all SNPs sharing the same database parent
    /// </summary>
    public class SibblingSNPs
    {
        private string databaseName;

        private List<SNPSite> all;

        public SibblingSNPs(List<SNPSite> sites)
        {
            if (sites.Count > 0) this.databaseName = sites[0].DatabaseParent;
            else this.databaseName = null;
            this.all = sites;
        }

        /// <summary>
        /// Get the number of tags having SNPs 
        /// </summary>
        public int Count_TagsHavingSNPs
        {
            get
            {
                Dictionary<string, int> dict = new Dictionary<string, int>();
                List<SNP> subsNPs = this.SNPs_Valid;
                foreach (SNP s in subsNPs)
                {
                    if (!dict.ContainsKey(s.QueryParent)) dict.Add(s.QueryParent, 1);
                }
                return dict.Count;
            }

        }


        /// <summary>
        /// Return all SNP sites mapping to this database parent
        /// </summary>
        public List<SNPSite> SNPSites
        {
            get
            {

                return this.all;
            }
        }



        /// <summary>
        /// Return all valid SNPs mapping to this database parent
        /// </summary>
        public List<SNP> SNPs_Valid
        {
            get
            {
                return SNPSiteUtility.GetSNPs_Valid(all);
            }
        }

        //Return all valid and invalid SNPs mapping to this database parent
        public List<SNP> SNPs_ValidAndInvalid
        {
            get
            {
                return SNPSiteUtility.GetSNPs_ValidAndInvalid(all);
            }
        }

        /// <summary>
        /// Get the number of transitions at this database parent
        /// </summary>
        public int Transitions
        {
            get
            {
                int trans = 0;
                foreach (SNP s in this.SNPs_Valid)
                {
                    if (s.Transition == true) trans++;
                }
                return trans;
            }
        }

        /// <summary>
        /// Get the number of transversions at this datbase parent
        /// </summary>
        public int Transversions
        {
            get
            {
                int count = 0;
                foreach (SNP s in this.SNPs_Valid)
                {
                    if (s.Transversion == true) count++;
                }
                return count;
            }
        }

        /// <summary>
        /// Get a subset of snps
        /// </summary>
        /// <param name="minimumCount"></param>
        /// <param name="minimumProminenceInPercent"></param>
        public SibblingSNPs GetSubset(int minimumCount, float minimumProminenceInPercent)
        {
            List<SNPSite> sites = new List<SNPSite>();

            foreach (SNPSite s in this.all)
            {
                if (s.Count_Valid >= minimumCount && s.SNPSiteProminence_Percent >= minimumProminenceInPercent) sites.Add(s);
            }
            return new SibblingSNPs(sites);
        }

        public SibblingSNPs GetSubset(bool showPlusPlus, bool showPlusMinus, bool show_A, bool show_T, bool show_C, bool show_G, bool showIndel, int minCountAtLoci, int minAlleleAtLoci, int maxAlleleAtLoci, double minPicAtLoci, int minimumDistanceFromAlignmentsEnd, int? maxLowQualityTokens, int? minSequenceQualitySite, int? minSequenceQualityNeighborhood)
        {
            //First update the SNPLoci in all, that is create a subset of the SNPs in all
            List<SNPSite> tempSites = new List<SNPSite>();
            foreach (SNPSite locus in all)
            {
                SNPSite s = locus.GetSubset(showPlusPlus, showPlusMinus, show_A, show_T, show_C, show_G, showIndel, minimumDistanceFromAlignmentsEnd, maxLowQualityTokens, minSequenceQualitySite, minSequenceQualityNeighborhood);
                if (s.Count_Valid > 0) tempSites.Add(s);
            }

            List<SNPSite> toBuild = new List<SNPSite>();
            foreach (SNPSite locus in tempSites) //then  check if the subset fullfills the other  requirements
            {
                if (locus.Count_Alleles >= minAlleleAtLoci && locus.Count_Alleles <= maxAlleleAtLoci && locus.Count_Valid >= minCountAtLoci && locus.PIC > (minPicAtLoci - 0.001)) toBuild.Add(locus);
            }
            return new SibblingSNPs(toBuild);
        }


        /// <summary>
        /// Get the name of the datbase parent
        /// </summary>
        public string DatabaseParent
        {
            get
            {
                return this.databaseName;
            }
        }


        /// <summary>
        /// Get the number of SNP sites mapping to this database parent
        /// </summary>
        public int Count_SNPSites
        {
            get
            {
                return this.all.Count;
            }
        }



        /// <summary>
        /// Get the number of SNPs mapping to this database parent
        /// </summary>
        public int Count_SNPs_Valid
        {
            get
            {
                return this.SNPs_Valid.Count;
            }

        }

        /// <summary>
        /// Get the number of SNPs being derived from the plus/plus alignment (sense)
        /// </summary>
        public int Count_PlusPlus
        {
            get
            {
                int count = 0;
                foreach (SNP s in this.SNPs_Valid)
                {
                    if (s.PlusPlusStrand == true) count++;
                }
                return count;
            }
        }





        public static string ToStringHeader(string separatingString)
        {
            return String.Format("gene_ID{0}count_SNPs_all{0}count_tags_having_a_SNP_all{0}count_SNP_sites_all{0}count_SNPs_subset{0}count_tags_having_a_SNP_subset{0}count_SNP_sites_subset{0}sense_percent{0}transitions_subset{0}transversions_subset", separatingString);
        }



        public string ToString_All(string separatingString)
        {
            return this.Count_SNPs_Valid.ToString() + separatingString
            + this.Count_TagsHavingSNPs.ToString() + separatingString
             + this.Count_SNPSites.ToString();
        }


        public string ToString_Subset(string separatingString)
        {
            StringBuilder sb = new StringBuilder("");

            sb.Append(this.Count_SNPs_Valid.ToString() + separatingString);
            sb.Append(this.Count_TagsHavingSNPs.ToString() + separatingString);
            sb.Append(this.Count_SNPSites.ToString() + separatingString);
            sb.Append(String.Format("{0:f}", 100.0 * (double)this.Count_PlusPlus / (double)this.Count_SNPs_Valid) + separatingString);
            sb.Append(this.Transitions.ToString() + separatingString);
            sb.Append(this.Transversions.ToString());
            return sb.ToString();
        }

        public string ToString(string separatingString)
        {


            return this.DatabaseParent + separatingString
            + ToString_All(separatingString) + separatingString
            + ToString_Subset(separatingString);
        }

    }





}



