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
using Bio.SNP;
using System.IO;

namespace Bio.SNP.IO
{

    /// <summary>
    /// Reads SSRs from a TabDelimitedSSR-File
    /// </summary>
    public class RefSNPReader
    {
        private StreamReader sr;
        public static int linecount = 0;

        public RefSNPReader(StreamReader sr)
        {
            this.sr = sr;

        }

        public RefSNP GetNextRefSNP()
        {

            string line = "";
            //As long as there are lines in rhe file proceed!
            while (((line = sr.ReadLine()) != null))
            {
                linecount++;
                //Get rid of space
                line = line.Trim();
                if (line == "" || linecount < 6) continue;
                else
                {
                    return ParseTabDelimitedRefSNP(line);

                }
            }

            return null;
        }


        public List<RefSNP> GetRefSNPs(int number)
        {
            List<RefSNP> toRet = new List<RefSNP>();
            RefSNP temp;
            while (toRet.Count<number && (temp = GetNextRefSNP()) != null)
            {
                toRet.Add(temp);
            }

            return toRet;
        }


        public List<RefSNP> GetAllRefSNPs()
        {
            List<RefSNP> toRet = new List<RefSNP>();
            RefSNP temp;
            while ((temp = GetNextRefSNP()) != null)
            {
                toRet.Add(temp);
            }

            return toRet;
        }





        private RefSNP ParseTabDelimitedRefSNP(string toParse)
        {
            /*
* Column   Data
1      RefSNP id (rs#)
2      mapweight where
    1 = Unmapped
    2 = Mapped to single position in genome
    3 = Mapped to 2 positions on a single chromosome
    4 = Mapped to 3-10 positions in genome (possible paralog hits)
    5 = Mapped to >10 positions in genome
3      snp_type where
    0 = Not withdrawn.
    1 = Withdrawn. There are several reasons for withdrawn, the
        withdrawn status is fully defined in the asn1, flatfile,
        and XML descriptions of the RefSNP. See /specs/docsum_2005.asn
        for a full definition of snp-type values.
4      Total number of chromosomes hit by this RefSNP during mapping
5      Total number of contigs hit by this RefSNP during mapping
6      Total number of hits to genome by this RefSNP during mapping
7      Chromosome for this hit to genome
8      Contig accession for this hit to genome
9      Version number of contig accession for this hit to genome
10      Contig ID for this hit to genome
11      Position of RefSNP in contig coordinates
12      Position of RefSNP in chromosome coordinates (used to order report)
    Locations are specified in NCBI sequence location convention where:
           x, a single number, indicates a feature at base position x
           x..y, denotes a feature that spans from x to y inclusive.
           x^y, denotes a feature that is inserted between bases x and y
13      Genes at this same position on the chromosome
14      Average heterozygosity of this RefSNP
15      Standard error of average heterozygosity
16      Maximum reported probability that RefSNP is real. (For computationally-
     predicted submissions)
17      Validated status
     0 = No validation information
     1 = Cluster has 2+ submissions, with 1+ submission assayed 
         with a non-computational method
     2 = At least one subsnp in cluster has frequency data submitted
     3 = Non-computational method in cluster and frequency data present
     4 = At lease one subsnp in cluster has been experimentally 
         validated by submitter
         for other validation status values, please see:
         <a href="ftp://ftp.ncbi.nih.gov/snp/database/organism_shared_data/SnpValidationCode.bcp.gz">ftp://ftp.ncbi.nih.gov/snp/database/organism_shared_data
         /SnpValidationCode.bcp.gz</a>
18      Genotypes available in dbSNP for this RefSNP
     1 = yes
     0 = no
19      Linkout available to submitter website for further data on the RefSNP
     1 = yes
     0 = no
20      dbSNP build ID when the refSNP was first created (i.e. the create date)
21      dbSNP build ID of the most recent change to the refSNP cluster. The
 date of the change is represented by the build ID which has an
 approximate date/time associated with it. (see:
 http://www.ncbi.nlm.nih.gov/projects/SNP/buildhistory.cgi)
22       Mapped to a reference or alternate (e.g. Celera) assembly 
*/
            string[] division = toParse.Split(new char[] { '\t' });
            if (division.Length != 22) throw new Exception("Error during reading of Tab Delimited RefSNP-File.");
            division[10] = division[10].Trim();
            division[11] = division[11].Trim();
            division[13] = division[13].Trim();
            division[14] = division[14].Trim();
            division[15] = division[15].Trim();
            System.Globalization.NumberFormatInfo ni = new System.Globalization.NumberFormatInfo();
            ni.NumberDecimalSeparator = ".";
#if DEBUG
            Convert.ToInt32(division[0]);
            Convert.ToInt32(division[1]); //int mapweight
            Convert.ToInt32(division[2]); //int snpType
            Convert.ToInt32(division[3]); //int chromosomeHits
            Convert.ToInt32(division[4]); //int contigHits
            Convert.ToInt32(division[5]); //int genomeHits
            string temp = division[6]; //string chromosome
            temp = division[7]; //string contigAccession
            Convert.ToInt32(division[8]); //int contigVersionNumber
            temp = division[9]; //string contigID
            int? testint = division[10] == "" ? (int?)null : Convert.ToInt32(division[10]);//int startPostionInContig
            testint = division[11] == "" ? (int?)null : Convert.ToInt32(division[11]);//int startPostionInChromosome
            string[] temp2 = division[12].Contains(",") ? division[12].Split(new char[] { ',' }) : null;//string[] genes
            double? av = division[13] == "" ? (double?)null : Double.Parse(division[13], ni); //double? averageHeterozygosity
            av = division[14] == "" ? (double?)null : Double.Parse(division[14], ni); //double? stdDev_averageHeterozygosity
            av = (division[15] == "") ? (double?)null : Double.Parse(division[15], ni); //double? p_realRefSnP
            Convert.ToInt32(division[16]);     //int validationStatus
            bool test = division[17] == "1" ? true : false;      //bool genotypesAvailable
            test = division[18] == "1" ? true : false;      //bool linkOutAvailable
            Convert.ToInt32(division[19]);      //int dbSNPBuildFirtsBuild
            Convert.ToInt32(division[20]);     //int dbSNPBuildRecentChange
            temp = division[21];
#endif



            //Add the Microsatellite
            return new RefSNP(Convert.ToInt32(division[0]) //int RefSNPId
                             , Convert.ToInt32(division[1]) //int mapweight
                             , Convert.ToInt32(division[2]) //int snpType
                             , Convert.ToInt32(division[3]) //int chromosomeHits
                             , Convert.ToInt32(division[4]) //int contigHits
                             , Convert.ToInt32(division[5]) //int genomeHits
                             , division[6]                  //string chromosome
                             , division[7]                  //string contigAccession
                             , Convert.ToInt32(division[8]) //int contigVersionNumber
                             , division[9]                  //string contigID
                             , division[10] == "" ? (int?)null : Convert.ToInt32(division[10])//int startPostionInContig
                             , division[11] == "" ? (int?)null : Convert.ToInt32(division[11])//int startPostionInChromosome
                             , division[12].Contains(",") ? division[12].Split(new char[] { ',' }) : null//string[] genes
                             , division[13] == "" ? (double?)null : Double.Parse(division[13], ni) //double? averageHeterozygosity
                             , division[14] == "" ? (double?)null : Double.Parse(division[14], ni) //double? stdDev_averageHeterozygosity
                             , division[15] == "" ? (double?)null : Double.Parse(division[15], ni) //double? p_realRefSnP
                             , Convert.ToInt32(division[16])     //int validationStatus
                             , division[17] == "1" ? true : false      //bool genotypesAvailable
                             , division[18] == "1" ? true : false      //bool linkOutAvailable
                             , Convert.ToInt32(division[19])      //int dbSNPBuildFirtsBuild
                             , Convert.ToInt32(division[20])     //int dbSNPBuildRecentChange
                             , division[21]                      //string mappedToReference
                             );

        }

    }


    /// <summary>
    /// Instance of a SNP filewriter; SNPs are written in the tab delimited format
    /// </summary>
    public class SNPWriter
    {
        private StreamWriter sw;
        private bool wroteIni = false;
        private bool wroteGeneralInfo = false;
        public SNPWriter(StreamWriter sw)
        {
            this.sw = sw;

        }
        public void WriteSNPs(List<SNP> snps)
        {

#if DEBUG
            System.Diagnostics.Stopwatch stopw = new System.Diagnostics.Stopwatch();
            stopw.Start();
#endif


            foreach (SNP s in snps)
            {

                WriteSNP(s);
            }

#if DEBUG
            stopw.Stop();
            System.Diagnostics.Debug.WriteLine(String.Format("{0} SNPWriter; {1} SNPs writen in TabDelimited style output format; Time used: {2} ms", DateTime.Now.TimeOfDay, snps.Count, stopw.ElapsedMilliseconds));
#endif
        }


        public void WriteSNP(SNP snp)
        {
            if (!wroteIni) //tdatabase_name\tquery_name\tstrand\tdatabase_char\tquery_char\tposition_database\tposition_query\tallele_identification\tinvalid_allele\tallele_id\tsynonymous")
            {
                sw.WriteLine("#");
                sw.WriteLine("#\tdatabase_name\tquery_name\tstrand\tdatabase_char\tquery_char\tsiteID_database\tsiteID_query\tvalid\tdistance_from_alignment_ends\talignment_Quality_Token\tsequence_quality_at_site\tsequence_quality_in_neighborhood\tallele_id");
                sw.WriteLine("#");
                sw.WriteLine();
                wroteIni = true;
            }

            if (!wroteGeneralInfo && snp !=null)
            {

                sw.WriteLine(String.Format(";Mode {0}", snp.IsFullSNP == true ? "fullSNP" : "halfSNP"));
                sw.WriteLine();
                wroteGeneralInfo = true;
            }

            sw.WriteLine(String.Format("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}"
    , snp.DatabaseParent.Replace('\t', ' ')
    , snp.QueryParent.Replace('\t', ' ')
    , snp.PlusPlusStrand == true ? "+/+" : "+/-"
    , snp.DatabaseChar
    , snp.QueryChar
    , snp.SiteID_Database
    , snp.SiteID_Query
    , snp.Valid.ToString()
    , snp.DistanceFromAlignmentEnd
    , snp.AlignmentQualityTokens
    , snp.SequenceQualityAtSite
    , snp.SequenceQualityInNeighborhood
    , snp.AlleleID != null ? snp.AlleleID.Replace('\t', ' ') : ""
    ));
        }
    }


    /// <summary>
    /// Reads SSRs from a TabDelimitedSSR-File
    /// </summary>
    public class SNPReader
    {
        private StreamReader sr;
        private bool? snpAlleleMode = null;


        public SNPReader(StreamReader sr)
        {
            this.sr = sr;

        }

        public SNP GetNextSNP()
        {
            string line = "";

            //As long as there are lines in rhe file proceed!
            while (((line = sr.ReadLine()) != null))
            {
                //Get rid of space
                line = line.Trim();
                if (line == "" || line.StartsWith("#")) continue;
                else if (snpAlleleMode == null && line.StartsWith(";"))
                {
                    //sw.WriteLine(String.Format(";Mode {0}",snps[0].IsSNPAllele==true?"SNP-Alleles":"SNPs"));
                    if (line.StartsWith(";Mode"))
                    {
                        string test = line.Split(new char[] { ' ' })[1];
                        if (test == "halfSNP") this.snpAlleleMode = false; //(String.Format(";Mode {0}",snps[0].IsFullSNP==true?"fullSNP":"halfSNP"));
                        else if (test == "fullSNP") this.snpAlleleMode = true;
                        else throw new InvalidDataException("Invalid SNP mode");
                    }
                }
                else
                {
                     return  ParseTabDelimitedSNP(line);
                }

            }


            return null;
        }



        public List<SNP> GetAllSNPs()
        {

            List<SNP> toRet = new List<SNP>();
            SNP temp;
            while ((temp = GetNextSNP()) != null)
            {
                toRet.Add(temp);
            }

            return toRet;

        }

        /// <summary>
        /// Retrieve a collection of SNPs which contains the specified number of SNPs
        /// </summary>
        /// <param name="number">the total number of SNPs to be returned</param>
        /// <returns></returns>
        public List<SNP> GetNextSNPs(int number)
        {

            List<SNP> toRet = new List<SNP>();
            SNP temp;
            while (toRet.Count<number && (temp = GetNextSNP()) != null)
            {
                toRet.Add(temp);
            }

            return toRet;
        }



        //Extracts the SSR related information from the sputnikInfo string and stores the SSR in the database
        private SNP ParseTabDelimitedSNP(string toParse)
        {
            //      0            1           2          3          4                    5               6        7                    8                      9                      10                      11                               12
            //database_name\tquery_name\tstrand\tdatabase_char\tquery_char\tposition_database\tposition_query\tvalid\tdistance_from_alignment_ends\talignment_Quality_Token\tsequence_quality_at_site\tsequence_quality_in_neighborhood\tallele_id");

            string[] division = toParse.Split(new char[] { '\t' });
            if (division.Length < 7) throw new Exception("Error during reading of Tab Delimited SNP-File.");


            SNP s = new SNP(division[0]
                , division[1]
                , division[3][0]
                , division[4][0]
                , division[5]
            , division[6]);

            s.PlusPlusStrand = division[2] == "+/+" ? true : false;
            s.IsFullSNP = this.snpAlleleMode.Value;
            s.Valid = division[7] == "True" ? true : false;


            if (division.Length > 8)
            {
                if (division[8] != null && division[8] != "") s.DistanceFromAlignmentEnd = Convert.ToInt32(division[8]);
            }

            if (division.Length > 9)
            {
                if (division[9] != null && division[9] != "") s.AlignmentQualityTokens = Convert.ToInt32(division[9]);
            }
            if (division.Length > 10)
            {
                if (division[10] != null && division[10] != "") s.SequenceQualityAtSite = Convert.ToInt32(division[10]);
            }
            if (division.Length > 11)
            {
                if (division[11] != null && division[11] != "") s.SequenceQualityInNeighborhood = Convert.ToInt32(division[11]);
            }
            if (division.Length > 12)
            {
                division[12] = division[12].Trim();
                if (division[12] != null && division[12] != "") if (division[12] != "") s.AlleleID = division[12];
            }

            //Add the Microsatellite
            return s;

        }

    }


    public class SNPSiteIO
    {
        /// <summary>
        /// Write the Loci into a file, info for the subset of SNPs at the locus
        /// </summary>
        /// <param name="sw"></param>
        /// <param name="snpSites"></param>
        public static void WriteSNPSites(StreamWriter sw, List<SNPSite> snpSites)
        {
            sw.WriteLine("#SNP LOCI");
            sw.WriteLine("#database_parent\tdatabase_position\tcount_subset\tcount_alleles_subset\tcount_pic_subset");
            sw.WriteLine("#only loci which having a subset count larger than zero are writen");
            sw.WriteLine();

            foreach (SNPSite s in snpSites)
            {
                if (s.Count_Valid > 0) //Write the SNP site only if the counts are larger than zero
                {
                    sw.WriteLine(String.Format("{0}\t{1}\t{2}\t{3}\t{4:f}"
                        , s.DatabaseParent
                        , s.SiteID_Database
                        , s.Count_Valid
                        , s.Count_Alleles
                        , s.PIC));
                }
            }
        }


        public static List<SNPSite> ReadSNPSites(StreamReader sr)
        {
            List<SNPSite> snpSites = new List<SNPSite>();
            string line = "";

            //As long as there are lines in rhe file proceed!
            while ((line = sr.ReadLine()) != null)
            {
                //Get rid of space
                line = line.Trim();
                if (line == "" || line.StartsWith("#")) continue;
                else
                {
                    //              0                   1               2           3                           4                   5           6               7           8
                    // ("#database_parent\tdatabase_position\tcount_subset\tcount_alleles_subset\tcount_pic_subset");
                    string[] split = line.Split(new char[] { '\t' });
                    SNPSite s = new SNPSite(split[0], split[1]);
                    snpSites.Add(s);
                }
            }

            return snpSites;

        }
    }


}