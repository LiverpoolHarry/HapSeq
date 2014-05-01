/*
 * The MIT License
 *
 * Copyright (c) 2009 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */
package src.java.net.sf.picard.analysis;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import src.java.net.sf.picard.cmdline.Option;
import src.java.net.sf.picard.io.IoUtil;
import src.java.net.sf.picard.reference.ReferenceSequence;
import src.java.net.sf.samtools.SAMFileHeader;
import java.io.File;

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.text.DateFormat;
import java.text.DecimalFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Date;
import java.util.HashMap;
import java.util.Iterator;
import java.util.TreeMap;
import src.java.net.sf.samtools.SAMRecord;

/**
 * Gets total length of contigs from BAM / SAM file
 * @author harry
 *
 */
public class ContigMetrics extends SinglePassSamProgram {

   private File tempFile = new File("temp.tmp");// default for file names so that they are not required
   @Option(doc = "Boolean. If set to true calculate mean quality over aligned reads only")
   public boolean ALIGNED_READS_ONLY = true;
   @Option(shortName = "PF", doc = "Boolean. If set to true calculate mean quality over PF reads only")
   public boolean PF_READS_ONLY = false;
   @Option(shortName = "MRL", doc = "Integer. Minimum allowed read length. Many reads may be less than the nominal read length. Default 50")
   public int MIN_READ_LENGTH = 10;
   @Option(shortName = "MPL", doc = "Integer. Maximum allowed length of Paired reads minimum position of forward read to maximum position of reverse read. Default 1000")
   public int MAX_PAIR_LENGTH = 1000;
   @Option(shortName = "MCL", doc = "Integer. Minimum allowed contig length other wise many contigs may be less than the nominal read length. Default 50")
   public int MIN_CONTIG_LENGTH = 50;
   @Option(shortName = "LG", doc = "Integer. Maximum length of gap between contigs to use for calculating mean and standard deviation of gaps. Large gaps would skew values. ")
   public int LENGTH_OF_GAP = 1000;
   @Option(shortName = "MQ", doc = "Integer. Minimum map quality to accept a read. Reads with less than this value will be excluded from analysis. Set to 0 to include all reads")
   public int MIN_MAP_QUAL = 20;
   @Option(shortName = "MZG", doc = "Double. Number of standard deviations from median gap between contigs to use for building scaffolds.")
   public Double MAX_Z_GAP = 3.0;
   @Option(shortName = "MAG", doc = "Integer. Maxiumum length of gap between contigs to use for building scaffolds.")
   public Integer MAX_ABS_GAP = 500;
   @Option(shortName = "PSC", doc = "Double. Proportion of total scaffold length that is expected to be covered by contigs. Default 1.0. With default no scaffolds are constructed by this strategy")
   public Double PROP_SCAFFOLD_COVERED = 1.0;
   @Option(shortName = "MS", doc = "Integer. Minimum scaffold length for inclusion in metric calculations and interval list out put")
   public Integer MIN_SCAFFOLD = 2000;
   @Option(shortName = "ILZ", doc = "String. Name for interval list output file of scaffold coordinates in format chr:start-end. Must have '.interval_list' extension if used for Single Indivdual Haplotyper. Scaffolds built with gaps between contigs less than mean + MAX_Z_GAP sigma")
   public File INTERVAL_LIST_ZSCORE = tempFile;
   @Option(shortName = "RF", doc = "String. Name of input file of repeat coordinates in format chr:start-end. If present will be used to exclude scaffolds that are within a single repeat.")
   public File REPEATS_FILE = tempFile;
   @Option(shortName = "ILP", doc = "String. Name for interval list output file of scaffold coordinates in format chr:start-end. Must have '.interval_list' extension if used for Single Indivdual Haplotyper. Scaffolds built with total gaps between contigs up to total length of target with zero coverage. ")
   public File INTERVAL_LIST_PROPORTION = tempFile;
   @Option(shortName = "CGH", doc = "String. Name of output file for histogram of gap lengths.")
   public File CONTIG_GAPS_FILE = tempFile;
   @Option(shortName = "CF", doc = "String. Name of output file for list of contigs ")
   public File CONTIGS_FILE = tempFile;
   @Option(shortName = "CLH", doc = "String. Name of output file for histogram of contig lengths. ")
   public File CONTIG_LENGTHS_FILE = tempFile;
   @Option(shortName = "COV", doc = "String. Name of output file for histogram of coverage. ")
   public File COVERAGE_FILE = tempFile;
   @Option(shortName = "CCIS", doc = "String. Name of output file for histogram of counts of contigs in scaffolds. ")
   public File COUNT_CONTIGS_IN_SCAFFOLDS_FILE = tempFile;
   private int totalContigLength = 0;//Cumulative length of contigs
   private long thisContig = 0;
   private long contigStart = 0;
   private ArrayList<Coordinates> readPairs;
   private HashMap<String, Integer> readPairStarts;//hold starts until paired end is found
   //private ArrayList<String> readNames;//list of read pair names in position order
   private int contigCount = 0;
   private int shortContigs = 0;
   private float meanShortContigLength = 0;
   private int contigReadCount = 0;
   private ArrayList<Coordinates> contigs;
   private String previousReferenceName = "XX";
   private String refName = "XX";
   private int countBelowMinMapQual = 0;
   private long previousEnd = 0; //End position of previous record
   private long previousStart = 0; //Start position of previous record
   private String previousRefName = "XX";
   private long readLength = 0; //rec.getReadLength() which seems to equal original read length
   private long totalReadLength = 0;//sum of rec.getReadLength() which seems to equal original read length
   private long totalObservedReadLength = 0;//sum of (end - start) which looks like mapped read length
   private long totalReadCount = 0;
   private ReferenceSequence refSeq;
   private String referenceName;
   private long scaffoldStart = 0;
   private long thisScaffold = 0;
   private int scaffoldContigCount = 0;
   private HashMap<Integer, Long> percentiles;
   private HashMap<Long, Integer> contigGaps;
   private HashMap<Long, Integer> contigLengths;
   private HashMap<Long, Integer> countContigsInScaffolds;
   private HashMap<Long, Integer> coverage; //Whole chromosome coverage; key=coverage; value = number of bases at key coverage
   private HashMap<Long, Long> contigCoverage; //Current contig coverage; key=coverage; value = number of bases at key coverage
   private HashMap<String, TreeMap<Long, Long>> chrMaps;// A treeMap of scaffold starts and ends for each chromosome;
   private long gapLength;//gaps less than LENGTH_OF_GAP
   private long allGapLength;//all gaps
   private int gapCount;//count gaps less than LENGTH_OF_GAP
   private ArrayList<Long> gaps;
   private ArrayList<Long> allGaps;
   private BufferedWriter writer;
   private BufferedWriter gap_length_writer;
   private String inputFilename;
   private DecimalFormat myFormatter;
   private DecimalFormat floatFormat;
   private DateFormat dateFormat;
   private int pairedReads = 0;
   private int firstOfPairCount = 0;
   private int singletonReads = 0;
   private Double fractionOfContigLensAboveN50;
   private long maxContigLength = 0; //length of longest observed contig
   private Double readFrac = 0.01;//added to position to make unique position allowing for multiple reads mapping to same position

   /** Required main method. */
   public static void main(final String[] args) {
      System.exit(new ContigMetrics().instanceMain(args));
   }

   @Override
   protected void setup(final SAMFileHeader header, final File samFile) {
      IoUtil.assertFileIsWritable(OUTPUT);
      writer = IoUtil.openFileForBufferedWriting(OUTPUT, true);

      SAMFileHeader outHeader = new SAMFileHeader();
      inputFilename = INPUT.getName();
      inputFilename = inputFilename.substring(0, inputFilename.length() - 4);
      //pad filename with spaces so that output always aligns
      inputFilename = String.format("%1$-20s", inputFilename);
      gaps = new ArrayList<Long>();
      allGaps = new ArrayList<Long>();
      contigs = new ArrayList<Coordinates>();
      myFormatter = new DecimalFormat("###,###");
      floatFormat = new DecimalFormat("###,###.###");
      dateFormat = new SimpleDateFormat("yyyy/MM/dd HH:mm:ss");
      percentiles = new HashMap<Integer, Long>(101);
      contigGaps = new HashMap<Long, Integer>();
      contigLengths = new HashMap<Long, Integer>();
      contigCoverage = new HashMap<Long, Long>();
      coverage = new HashMap<Long, Integer>();
      readPairs = new ArrayList<Coordinates>();
      readPairStarts = new HashMap<String, Integer>();
      countContigsInScaffolds = new HashMap<Long, Integer>();
      chrMaps = new HashMap<String, TreeMap<Long, Long>>();
   }

   @Override
   protected void acceptRead(final SAMRecord rec, final ReferenceSequence ref) {
      // Skip unwanted records
      refSeq = ref;
      if (PF_READS_ONLY && rec.getReadFailsVendorQualityCheckFlag()) {
         return;
      }
      if (ALIGNED_READS_ONLY && rec.getReadUnmappedFlag()) {
         return;
      }
      if (rec.getNotPrimaryAlignmentFlag()) {
         return;
      }
      referenceName = rec.getReferenceName();
      if (previousRefName.contentEquals("XX")) {
         previousRefName = referenceName;
      }
      //In case running on sample with multiple chromsosomes but not been tested like this
      if (!referenceName.contentEquals(previousRefName)) {
         buildContigs();
         Date date = new Date();
         System.err.println(dateFormat.format(date) + "  " + inputFilename + " Processed " + previousRefName);
         previousRefName = rec.getReferenceName();
      }

      totalReadCount++;
      if (rec.getMappingQuality() < MIN_MAP_QUAL) {
         countBelowMinMapQual++;
         return;
      }
      String readName = rec.getReadName();
      //create uniqe numerical index by adding a fractional increment to start position
      Double index = rec.getAlignmentStart() + 0.0;
      if (rec.getAlignmentStart() == previousStart) {
         index = previousStart + readFrac;
         readFrac += 0.01;
      }
      else {
         readFrac = 0.01;
      }

      if (rec.getReadPairedFlag()) {
         if (rec.getFirstOfPairFlag()) {
            //check if pairs are on same chromosome and discard read if they are not
            if (!rec.getMateReferenceName().contentEquals(rec.getReferenceName())) {
               return;
            }
            firstOfPairCount++;
         }

         pairedReads++;
         if (rec.getFirstOfPairFlag()) {
            readPairStarts.put(readName, rec.getAlignmentStart());
            return;
         }

         if (rec.getSecondOfPairFlag()) {
            if (readPairStarts.containsKey(readName)) {
               int pairLength = rec.getAlignmentEnd() - readPairStarts.get(readName);
               if (pairLength >= MIN_READ_LENGTH && pairLength < MAX_PAIR_LENGTH) {
                  readPairs.add(new Coordinates(readName, referenceName,
                        readPairStarts.get(readName), rec.getAlignmentEnd(), index));
               }
            }
            return;
         }
         return;
      }
      else {
         int readSpan = rec.getAlignmentEnd() - rec.getAlignmentStart();
         if (readSpan >= MIN_READ_LENGTH && readSpan <= MAX_PAIR_LENGTH) {
            readPairs.add(new Coordinates(readName, referenceName,
                  rec.getAlignmentStart(), rec.getAlignmentEnd(), index));
            singletonReads++;
         }
      }
   }

   private void buildContigs() {
      Collections.sort(readPairs);
      for (Coordinates readPair : readPairs) {
         long start = readPair.getStart();
         long end = readPair.getEnd();
         refName = previousRefName;

         //swap start and end if end less than start
         if (end < start) {
            long tmp = start;
            start = end;
            end = tmp;
         }
         readLength = end - start + 1;
         totalReadLength += readLength;
         totalObservedReadLength += readLength;

         //set up initial values at start of run
         if (previousEnd == 0) {
            scaffoldStart = start;
            contigStart = start;
            previousEnd = end;
         }
         //assumes reads sorted
         //reads may have been trimmed and have a start after a previous read and end before a previous read
         //includes check for change of chromosome
         if (previousEnd < start) {
            thisContig = (previousEnd - contigStart + 1);
            //build hash of counts of Contig Lengths
            updateContigLengths();
            if (thisContig < MIN_CONTIG_LENGTH) {
               shortContigs++;
               meanShortContigLength += thisContig;
               previousEnd = end;
               contigStart = start;
            }
            else {
               contigCount++;
               totalContigLength += thisContig;
               Coordinates coords = new Coordinates("Contig " + contigCount, refName, contigStart, previousEnd);
               contigs.add(coords);

               //if contig still within scaffold
               Long gL = (start - previousEnd);
               allGaps.add(gL);
               allGapLength += gL;
               //build hash of counts by GAP lengths
               contigGaps.put(gL, (contigGaps.containsKey(gL) == false) ? 1 : contigGaps.get(gL) + 1);

               if (gL <= LENGTH_OF_GAP) {
                  gapLength += gL;
                  gaps.add(gL);
                  gapCount++;
                  scaffoldContigCount++;
               }
               //Add coverage of this contig to global stats
               updateCoverage(coords);
               contigReadCount = 1;
               contigStart = start;
            }
         }
         else {
            contigReadCount++;
         }
         setContigCoverage(start, end);

         //trimmed reads can mean the end might be before previous end
         if (end > previousEnd) {
            previousEnd = end;
         }

         //}
      }
      //Clear read pairs to save memory and prepare for next chromosome
      readPairs.clear();
   }

   private void updateCoverage(Coordinates coords) {
      ArrayList<Long> values = new ArrayList<Long>();
      //contigCoverage is hashmap of position=>coverage
      values.addAll(contigCoverage.values());
      int totalCoverage = 0;
      for (Long i : values) {
         Boolean keyExists = coverage.containsKey(i);
         coverage.put(i, (keyExists == false) ? 1 : coverage.get(i) + 1);
         totalCoverage += i;
      }
      Double meanCoverage = (1.0 * totalCoverage) / values.size();
      coords.setCoverage(meanCoverage);
      contigCoverage.clear();
   }

   protected void updateContigLengths() {
      contigLengths.put(thisContig, (contigLengths.get(thisContig) == null) ? 1
            : contigLengths.get(thisContig) + 1);
      if (maxContigLength < thisContig) {
         maxContigLength = thisContig;
      }
   }

   protected void setContigCoverage(long start, long end) {
      for (long i = start - contigStart; i <= end - contigStart; i++) {
         Boolean keyExists = contigCoverage.containsKey(i);
         contigCoverage.put(i, (keyExists == false) ? 1 : contigCoverage.get(i) + 1);
      }

   }

   //Get scaffolds based on a maximum gap length between contigs of mean + MAX_Z_GAP standard deviations
   private ScaffoldMetrics getScaffolds(long mgl) {
      Double stdDev = getStdDevgapLength(mgl);
      Double maxGap = MAX_Z_GAP * stdDev + mgl;
      ScaffoldMetrics sm = buildscaffolds(maxGap);
      sm.setStdev(stdDev);
      return sm;
   }

   //Get scaffolds based on a maximum absolute gap length between contigs
   private ScaffoldMetrics getAbsGapScaffolds() {
      Double stdDev = getStdDevgapLength(MAX_ABS_GAP);
      ScaffoldMetrics sm = buildscaffolds(MAX_ABS_GAP * 1.0);
      sm.setStdev(stdDev);
      return sm;
   }

   //Get scaffolds based on length of sequence with no coverage. The gaps are ordered and
   //then summed and the minimum gap length such that all gaps shorter than that add up to
   //the length of genome with no sequence is obtained.
   private ScaffoldMetrics buildScaffoldsBasedonProp() {
      //Calculate length of target present but with zero coverage
      Double totalAcceptedGapLength = (totalContigLength / PROP_SCAFFOLD_COVERED) - totalContigLength;
      Collections.sort(gaps);
      Double maxGap = 0.0; //Maximum gap length for joining contigs in scaffold
      Integer sumGaps = 0;
      Iterator i = gaps.iterator();
      while (sumGaps < totalAcceptedGapLength && i.hasNext()) {
         Integer g = (Integer) i.next();
         sumGaps += g;
         maxGap = g * 1.0;
      }
      ScaffoldMetrics sm = buildscaffolds(maxGap);
      return sm;
   }

   private ScaffoldMetrics buildscaffolds(Double maxGap) {
      long preEnd = contigs.get(0).getEnd();
      long scaStart = contigs.get(0).getStart();
      String preChr = contigs.get(0).getChr();
      long contCount = 0;
      int scaffCount = 0;
      long totalScaffLength = 0;
      ArrayList<Long> scaffLengths = new ArrayList<Long>();

      int allScaffCount = 0;
      int totalAllScaffLength = 0;
      ArrayList<Coordinates> scaffolds = new ArrayList<Coordinates>();
      ArrayList<Long> allScaffLengths = new ArrayList<Long>();
      Double scaffCoverage = 0.0;
      for (Coordinates contig : contigs) {
         if ((contig.getStart() - preEnd) > maxGap || !contig.getChr().contentEquals(preChr)) {
            if ((preEnd - scaStart) > MIN_SCAFFOLD) {
               Coordinates scaffold = new Coordinates("Scaffold " + scaffCount, preChr, scaStart, preEnd);
               scaffolds.add(scaffold);
               scaffold.setCoverage(scaffCoverage);
               Integer freq = countContigsInScaffolds.get(contCount);
               countContigsInScaffolds.put(contCount, (freq == null) ? 1 : freq + 1);
               scaffCoverage = 0.0;
            }
            scaStart = contig.getStart();
            scaffCoverage = contig.getCoverage();
         }
         else {
            contCount++;
            scaffCoverage += contig.getCoverage();
         }
         preEnd = contig.getEnd();
         preChr = contigs.get(0).getChr();
      }

      //Add final scaffold
      Coordinates scaffold = new Coordinates("Scaffold " + scaffCount, preChr, scaStart, preEnd);
      scaffolds.add(scaffold);

      ScaffoldMetrics sm = new ScaffoldMetrics(scaffolds, maxGap);//contCount, scaffCount, totalScaffLength, scaffLengths, allScaffCount, totalAllScaffLength, allScaffLengths, maxGap);
      sm.setContigCount(contCount);
      return sm;
   }

   private ScaffoldMetrics removeScaffoldsInRepeats(ScaffoldMetrics scaffs) {
      buildRepeatHash();
      ArrayList<Coordinates> filteredScaffs = new ArrayList<Coordinates>();
      for (Coordinates coords : scaffs.getScaffolds()) {
         //In a TreeMap floorKey is the highest key <= value
         if (chrMaps.get(coords.chr) != null) {
            if (chrMaps.get(coords.chr).floorKey(coords.start) != null) {
               Long maxRepeatStart = chrMaps.get(coords.chr).floorKey(coords.start);
               if (chrMaps.get(coords.chr).get(maxRepeatStart) < coords.end) {
                  filteredScaffs.add(coords);
               }
            }
         }
      }
      ScaffoldMetrics filtScaffs = new ScaffoldMetrics(filteredScaffs, scaffs.getMaxGap());

      if (filtScaffs.getAllScaffCount() > 0) {
         Long mean = (Long) (filtScaffs.getTotalAllScaffLength() / filtScaffs.getAllScaffCount());
         filtScaffs.setStdev(getStdDevgapLength(mean));
      }
      else {
         filtScaffs.setStdev(-1.0);
      }
      return filtScaffs;
   }

   //Create a seperate TreeMap for repeats in each chromosome and enter the TreeMap into
   // a HashMap (chrMaps) of TreeMaps
   // input file assumed to be in bed format
   protected void buildRepeatHash() {
      //Scaffold start end pairs
      TreeMap<Long, Long> repeatEnds = new TreeMap<Long, Long>();

      String preChr = "XYZ";

      try {
         BufferedReader in = new BufferedReader(new FileReader(REPEATS_FILE));
         String str;
         String chr = "Chr";
         int repeatCount = 0;
         while ((str = in.readLine()) != null) {
            String temp[] = str.split("\\s");
            chr = temp[0];
            Long start = Long.parseLong(temp[1]);
            Long end = Long.parseLong(temp[2]);
            if (!chr.contentEquals(preChr) && !preChr.contentEquals("XYZ")) {
               chrMaps.put(chr, (TreeMap<Long, Long>) repeatEnds.clone());
               repeatEnds.clear();
               repeatEnds.put(start, end);
            }
            else {
               repeatEnds.put(start, end);
            }
            preChr = chr;
         }
         chrMaps.put(chr, (TreeMap<Long, Long>) repeatEnds.clone());

         in.close();

      }
      catch (IOException e) {
         System.out.println("File not found: " + REPEATS_FILE);
      }


   }

   protected Double getStdDevgapLength(long mgl) {
      Double sumOfDiffs = 0.0;
      for (Long gap : gaps) {
         sumOfDiffs += (gap - mgl) * (gap - mgl);
      }

      Double stdDev = Math.sqrt(sumOfDiffs / gaps.size());
      return stdDev;

   }

   protected Long getN50contig() {
      ArrayList<Long> contLen = new ArrayList<Long>(contigLengths.keySet());
      long[] contLens = new long[contLen.size()];
      for (int i = 0; i < contLen.size(); i++) {
         contLens[i] = contLen.get(i);
      }
      Arrays.sort(contLens);
      int largest = contigLengths.get(contLens[0]) + 1;
      long runningTot = 0;
      int i = contLens.length - 1;
      long mid = (totalContigLength / 2);
      while (runningTot < mid && i > 0) {
         runningTot += contLens[i] * contigLengths.get(contLens[i]);
         largest = contigLengths.get(contLens[i]) + 1;
         i--;
      }
      //get exact N50
      Long diff = runningTot - mid;
      long lower = runningTot - contLens[i + 1] * contigLengths.get(contLens[i + 1]);
      long upper = runningTot;
      Long N50 = 0L;
      if (upper - lower == 0) {
         N50 = contLens[i];
      }
      else {
         N50 = contLens[i] + ((diff) / (upper - lower)) * (contLens[i + 1] - contLens[i]);
      }
      fractionOfContigLensAboveN50 = 1.0 * (i / contLens.length);
      return N50;
   }

   protected Long getMedian(HashMap<Long, Integer> hm) {
      Long median = 0L;
      ArrayList<Long> keys = new ArrayList<Long>(hm.keySet());
      Collections.sort(keys);
      Double midCount = 0.0;
      for (Long i : keys) {
         midCount += hm.get(i);
      }
      midCount = midCount / 2;
      Long runningTot = 0L;
      for (Long i : keys) {
         runningTot += hm.get(i);
         if (runningTot >= midCount) {
            median = i;
            break;
         }
      }
      return median;
   }

   protected String getPercentiles(ArrayList<Long> array, long totalLength) {
      Collections.sort(array);
      long runningTotal = 0;
      int percentile = 0;
      long inter = 0;
      int intervalCount = 0;
      int flag = 0;
      percentiles.clear();
      for (Long interval : array) {
         //need to multiply by 100.0 to ensure double is returned
         int currentPercentile = (int) ((100.0 * runningTotal) / totalLength);

         if (currentPercentile > percentile) {
            percentiles.put(percentile + 1, interval);
            printPercentile(percentile + 1, interval, runningTotal, intervalCount);
            intervalCount = 1;
            while (currentPercentile > percentile) {
               percentile++;
            }
         }
         runningTotal += interval;
         inter = interval;
         intervalCount++;

      }
      printPercentile(percentile, inter, runningTotal, intervalCount);
      return getN50(percentiles);
   }

   protected HashMap<String, Double> getCoverageStats(HashMap<Long, Integer> hash) {

      ArrayList<Long> keys = new ArrayList<Long>();
      keys.addAll(hash.keySet());
      // and sorting it
      Collections.sort(keys);
      Long runningTot = (long) 0;
      Long previousTot = (long) 0;
      Long previousKey = (long) 0;
      Long[] previousKeys = new Long[2];

      //Get N50 coverage
      for (Long key : keys) {
         if (runningTot < (totalReadLength / 2)) {
            previousKeys[0] = previousKey;
            previousKeys[1] = key;
            previousTot = runningTot;
            runningTot += key * hash.get(key);
            previousKey = key;
         }
      }
      long excess = runningTot - (totalReadLength / 2);
      Double fraction = (1.0) * excess / (runningTot - previousTot);
      Double n50Cov = 0.0;
      if (previousKeys[0] != null) {
         n50Cov = previousKeys[1] - (1 - (fraction * (previousKeys[1] - previousKeys[0])));
      }
      //Get standard deviation
      Double mean = 1.0 * totalReadLength / totalContigLength;
      Double stdDev = 0.0;
      for (Long k : keys) {
         stdDev += (k - mean) * (k - mean) * hash.get(k);
      }

      stdDev = java.lang.Math.sqrt(stdDev / totalContigLength);

      HashMap<String, Double> results = new HashMap<String, Double>(3);
      results.put("N50", n50Cov);
      results.put("stdDev", stdDev);

      return results;
   }

   @Override
   protected void finish() {
      buildContigs();
      System.err.println(inputFilename + " Processed " + previousRefName);

      //Tidy up last contig
      updateContigLengths();
      Coordinates coords = new Coordinates("Contig " + contigCount, referenceName, contigStart, previousEnd);
      contigs.add(coords);
      updateCoverage(coords);
      totalContigLength += (previousEnd - contigStart) + 1;
      contigCount++;
      meanShortContigLength = meanShortContigLength / shortContigs;

      int refSeqLength = refSeq.length();

      String meanCoverage = floatFormat.format(((totalReadLength * 1.0) / totalContigLength));
      String propRefCovered = floatFormat.format(((1.0 * totalContigLength) / refSeqLength));
      double meanScaffoldLength = 0.0;
      long meanGapLength = 0;
      long n50Gap = 0;
      HashMap<String, Double> coverageStats = getCoverageStats(coverage);


      System.err.println(inputFilename + " Read Count = " + myFormatter.format(totalReadCount));
      System.err.println(inputFilename + " Paired read count = " + myFormatter.format(pairedReads));
      System.err.println(inputFilename + " Singleton read count = " + myFormatter.format(singletonReads));
      System.err.println(inputFilename + " First read of pair count = " + myFormatter.format(firstOfPairCount));

      System.err.println(inputFilename + " Count reads below min map qual = " + myFormatter.format(countBelowMinMapQual));
      System.err.println(inputFilename + " Total reads Length = " + myFormatter.format(totalReadLength));
      System.err.println(inputFilename + " Total Contig Count = " + myFormatter.format(contigCount));
      System.err.println(inputFilename + " Total Contig Length = " + myFormatter.format(totalContigLength));
      System.err.println(inputFilename + " Length of longest Contig = " + myFormatter.format(maxContigLength));
      System.err.println(inputFilename + " Count contigs < " + MIN_CONTIG_LENGTH + " = " + myFormatter.format(shortContigs));
      System.err.println(inputFilename + " Mean length of contigs < " + MIN_CONTIG_LENGTH + " = " + myFormatter.format(meanShortContigLength));
      System.err.println(inputFilename + " N50 contig length = " + myFormatter.format(getN50contig()));
      System.err.println(inputFilename + " Median contig length = " + myFormatter.format(getMedian(contigLengths)));

      System.err.println(inputFilename + " Reference Length = " + myFormatter.format(refSeqLength));
      System.err.println(inputFilename + " Prop reference covered = " + propRefCovered);
      System.err.println(inputFilename + " Mean coverage = " + meanCoverage);
      System.err.println(inputFilename + " Standard deviation of mean Coverage = " + floatFormat.format(coverageStats.get("stdDev")));
      System.err.println(inputFilename + " N50 Coverage = " + floatFormat.format(coverageStats.get("N50")));

      String outString = "File=" + inputFilename + ";\nLENGTH_OF_GAP=" + LENGTH_OF_GAP +
            ";\nTotalContigLength=" + totalContigLength + ";\nContigCount=" + contigCount +
            ";\nTotalReadsLength=" + totalReadLength + ";\nTotalObservedReadsLength=" + totalObservedReadLength +
            ";\nReadCount=" + totalReadCount + ";\nReferenceLength=" + refSeqLength;

      if (gapCount > 0) {
         n50Gap = gaps.get((int) Math.floor(gaps.size() / 2));
         meanGapLength = (gapLength / gapCount);
         ScaffoldMetrics sZM = getScaffolds(meanGapLength);//scaffold metrics based on Z scores
         ScaffoldMetrics sAM = getAbsGapScaffolds();//Scaffold metrics based on absolute gap

         printHeader("Percentiles of Gap lengths");
         printHeader("Percentile\tInput File\tPercentile\tGap Length\tCumulative Gap Lengths\tCount");
         String n50gaps = getPercentiles(allGaps, allGapLength);

         System.err.println(inputFilename + " Mean Gap Length under " + LENGTH_OF_GAP + " = " + meanGapLength);
         System.err.println(inputFilename + " Standard deviation of gap Length under " + LENGTH_OF_GAP + " = " + floatFormat.format(sZM.getStdev()));
         System.err.println(inputFilename + " N50 Gap Length under " + LENGTH_OF_GAP + " = " + n50Gap);
         System.err.println(inputFilename + " N50 Gap Length: " + n50gaps);
         System.err.println(inputFilename + " Median Gap Length: " + myFormatter.format(getMedian(contigGaps)));
         System.err.println(inputFilename + " Scaffold Count = " + sAM.getScaffCount());

         if (sZM.getScaffCount() > 0) {
            meanScaffoldLength = (sZM.getTotalScaffLength() / sZM.getScaffCount());
            printHeader("Percentiles of Scaffold lengths");
            printHeader("Percentile\tInput File\tPercentile\tScaffold Length\tCumulative Scaffold Lengths\tCount");
            String n50scaffold = getPercentiles(sZM.getScaffLengths(), sZM.getTotalScaffLength());
            String n50AllScaffold = getPercentiles(sZM.getAllScaffLengths(), sZM.getTotalScaffLength());
            System.err.println(inputFilename + " Max gap between contigs to build scaffold mean + " + MAX_Z_GAP + " sigma = " + myFormatter.format(sZM.getMaxGap()));
            System.err.println(inputFilename + " Total Scaffold Length for scaffolds > " + MIN_SCAFFOLD + " from Z score = " + myFormatter.format(sZM.getTotalScaffLength()));
            System.err.println(inputFilename + " N50 Scaffold Length for scaffolds > " + MIN_SCAFFOLD + " from Z score " + n50scaffold);
            System.err.println(inputFilename + " Total Scaffold Length using gaps < mean + " + MAX_Z_GAP + " sigma = " + myFormatter.format(sZM.getTotalAllScaffLength()));
            System.err.println(inputFilename + " N50 Scaffold Length for all scaffolds from Z score " + n50AllScaffold);


            outString += ";\nTotal Gap Length Less Than LENGTH_OF_GAP=" + myFormatter.format(gapLength);
            outString += ";\nGap Count Less Than LENGTH_OF_GAP=" + myFormatter.format(gapCount);
            outString += ";\nMean Gap Length=" + myFormatter.format(meanGapLength);
            outString += ";\nTotal Scaffold Length=" + myFormatter.format(sZM.getTotalScaffLength());
            outString += ";\nScaffold Count Z score=" + myFormatter.format(sZM.getScaffCount());
            outString += ";\nMean Scaffold Length Z score=" + myFormatter.format(meanScaffoldLength);
            outString += ";\nN50 Scaffold Length for scaffolds > " + MIN_SCAFFOLD + " from Z score " + n50scaffold;
            outString += ";\nTotal Scaffold Length using gaps < mean + " + MAX_Z_GAP + " sigma = " + myFormatter.format(sZM.getTotalAllScaffLength());


            if (!INTERVAL_LIST_ZSCORE.getName().contentEquals("temp.tmp")) {
               printHistogram(countContigsInScaffolds, COUNT_CONTIGS_IN_SCAFFOLDS_FILE);
            }

            if (!REPEATS_FILE.getName().contentEquals("temp.tmp")) {
               ScaffoldMetrics filteredScaffs = removeScaffoldsInRepeats(sAM);


               if (filteredScaffs.getScaffCount() > 0) {
                  meanScaffoldLength = (filteredScaffs.getTotalScaffLength() / filteredScaffs.getScaffCount());
               }
               else {
                  meanScaffoldLength = 0.0;
               }
               printHeader("Percentiles of Scaffold lengths after scaffolds within repeats are filtered out");
               printHeader("Percentile\tInput File\tPercentile\tScaffold Length\tCumulative Scaffold Lengths\tCount");
               n50scaffold = getPercentiles(filteredScaffs.getScaffLengths(), filteredScaffs.getTotalScaffLength());
               n50AllScaffold = getPercentiles(filteredScaffs.getAllScaffLengths(), filteredScaffs.getTotalScaffLength());

               System.err.println(inputFilename + " Results after scaffolds within repeats have been filtered out");
               System.err.println(inputFilename + " Scaffold Count for filtered scaffolds > " + MIN_SCAFFOLD + " from with gaps < " + MAX_ABS_GAP + " = " + myFormatter.format(filteredScaffs.getScaffCount()));
               System.err.println(inputFilename + " Total Scaffold Length for filtered scaffolds > " + MIN_SCAFFOLD + " from with gaps < " + MAX_ABS_GAP + " = " + myFormatter.format(filteredScaffs.getTotalScaffLength()));
               System.err.println(inputFilename + " N50 Scaffold Length for filtered scaffolds > " + MIN_SCAFFOLD + " with gaps < " + MAX_ABS_GAP + " " + n50scaffold);
               System.err.println(inputFilename + " Scaffold Count  for filtered scaffolds with gaps < " + MAX_ABS_GAP + " = " + myFormatter.format(filteredScaffs.getAllScaffCount()));
               System.err.println(inputFilename + " Total Scaffold Length for filtered scaffolds with gaps < " + MAX_ABS_GAP + " = " + myFormatter.format(filteredScaffs.getTotalAllScaffLength()));
               System.err.println(inputFilename + " N50 Scaffold Length for all filtered scaffolds with gaps < " + MAX_ABS_GAP + " " + n50AllScaffold);

               outString += ";\nResults after scaffolds within repeats have been filtered out";
               outString += ";\nScaffold Count for filtered scaffolds > " + MIN_SCAFFOLD + " from with gaps < " + MAX_ABS_GAP + " = " + myFormatter.format(filteredScaffs.getScaffCount());
               outString += ";\nTotal Scaffold Length for filtered scaffolds > " + MIN_SCAFFOLD + " from with gaps < " + MAX_ABS_GAP + " = " + myFormatter.format(filteredScaffs.getTotalScaffLength());
               outString += ";\nN50 Scaffold Length for filtered scaffolds > " + MIN_SCAFFOLD + " with gaps < " + MAX_ABS_GAP + " " + n50scaffold;
               outString += ";\nScaffold Count  for filtered scaffolds with gaps < " + MAX_ABS_GAP + " = " + myFormatter.format(filteredScaffs.getAllScaffCount());
               outString += ";\nTotal Scaffold Length for filtered scaffolds with gaps < " + MAX_ABS_GAP + " = " + myFormatter.format(filteredScaffs.getTotalAllScaffLength());
               outString += ";\nN50 Scaffold Length for all filtered scaffolds with gaps < " + MAX_ABS_GAP + " " + n50AllScaffold;

               //remove white space from filename
               String samplename = inputFilename.trim();
               File filteredScaffoldsFile = new File(samplename + "_Gap_" + MAX_ABS_GAP + "_Filtered_Scaffolds.gff3");
               printCoordinates(filteredScaffoldsFile, FileType.gff3, filteredScaffs.getScaffolds());

            }
         }
         else {
            System.err.println("No scaffolds found");
         }

         if (sAM.getScaffCount() > 0) {

            String n50scaffold = getPercentiles(sAM.getScaffLengths(), sAM.getTotalScaffLength());
            String n50AllScaffold = getPercentiles(sAM.getAllScaffLengths(), sAM.getTotalScaffLength());

            outString += ";\nTotal Scaffold Length for scaffolds > " + MIN_SCAFFOLD + " and " + MAX_ABS_GAP + " absolute gap = " + myFormatter.format(sAM.getTotalScaffLength());
            outString += ";\nScaffold Count " + MAX_ABS_GAP + " absolute gap=" + myFormatter.format(sAM.getScaffCount());
            outString += ";\nTotal Scaffold Length using absolute gap <  " + MAX_ABS_GAP + "= " + myFormatter.format(sAM.getTotalAllScaffLength());
            outString += ";\nN50 Scaffold Length for scaffolds > " + MIN_SCAFFOLD + " from " + MAX_ABS_GAP + " absolute gap " + n50scaffold;
            outString += ";\nN50 Scaffold Length using absolute gap < " + MAX_ABS_GAP + " " + n50AllScaffold;

            System.err.println(inputFilename + " Total Scaffold Length " + MAX_ABS_GAP + " absolute gap=" + myFormatter.format(sAM.getTotalScaffLength()));
            System.err.println(inputFilename + " Scaffold Count " + MAX_ABS_GAP + " absolute gap=" + myFormatter.format(sAM.getScaffCount()));
            System.err.println(inputFilename + " Total Scaffold Length using absolute gap <  " + MAX_ABS_GAP + "= " + myFormatter.format(sAM.getTotalAllScaffLength()));
            System.err.println(inputFilename + " N50 Scaffold Length for scaffolds > " + MIN_SCAFFOLD + " with " + MAX_ABS_GAP + " absolute gap " + n50scaffold);
            System.err.println(inputFilename + " N50 Scaffold Length using absolute gap < " + MAX_ABS_GAP + " " + n50AllScaffold);

         }
         else {
            outString += ";\nNo scaffolds found using absolute gap " + MAX_ABS_GAP;
            System.err.println("No scaffolds found using absolute gap " + MAX_ABS_GAP);
         }

         if (PROP_SCAFFOLD_COVERED != 1.0) {
            ScaffoldMetrics pSC = buildScaffoldsBasedonProp();

            if (pSC.getScaffCount() > 0) {

               meanScaffoldLength = (pSC.getTotalScaffLength() / pSC.getScaffCount());
               printHeader("Percentiles of Scaffold lengths based on given proportion of target sequence with coverage (" + PROP_SCAFFOLD_COVERED + ")");
               printHeader("Percentile\tInput File\tPercentile\tScaffold Length\tCumulative Scaffold Lengths\tCount");
               String n50scaffold = getPercentiles(pSC.getScaffLengths(), pSC.getTotalScaffLength());
               String n50AllScaffold = getPercentiles(pSC.getAllScaffLengths(), pSC.getTotalScaffLength());
               System.err.println(inputFilename + " Max gap between contigs to build scaffold from proportion of coverage = " + myFormatter.format(pSC.getMaxGap()));
               System.err.println(inputFilename + " Total Scaffold Length for scaffolds > " + MIN_SCAFFOLD + " from proportion of coverage = " + myFormatter.format(pSC.getTotalScaffLength()));
               System.err.println(inputFilename + " Count Scaffolds > " + MIN_SCAFFOLD + " from proportion of coverage = " + myFormatter.format(pSC.getAllScaffCount()));
               System.err.println(inputFilename + " N50 Scaffold Length for scaffolds > " + MIN_SCAFFOLD + " from proportion of coverage = " + n50scaffold);
               System.err.println(inputFilename + " Total Scaffold Length from proportion of coverage = " + myFormatter.format(pSC.getTotalAllScaffLength()));
               System.err.println(inputFilename + " N50 Scaffold Length for all scaffolds from proportion of coverage = " + n50AllScaffold);

               if (!INTERVAL_LIST_PROPORTION.getName().contentEquals("temp.tmp")) {
                  printCoordinates(INTERVAL_LIST_PROPORTION, FileType.bed, pSC.getScaffolds());
               }

               outString =
                     outString + ";\nTotal Scaffold Length based on prop of target covered = " + myFormatter.format(pSC.getTotalScaffLength()) +
                     ";\nScaffold Count based on prop of target covered = " + myFormatter.format(pSC.getScaffCount()) +
                     ";\nMean Scaffold Length based on prop of target covered =" + myFormatter.format(meanScaffoldLength) +
                     ";\nMax gap between contigs to build scaffold from proportion of coverage = " + myFormatter.format(pSC.getMaxGap()) +
                     ";\nTotal Scaffold Length for scaffolds > " + MIN_SCAFFOLD + " from proportion of coverage = " + myFormatter.format(pSC.getTotalScaffLength()) +
                     ";\nCount Scaffolds > from proportion of coverage = " + myFormatter.format(pSC.getAllScaffCount()) +
                     ";\nN50 Scaffold Length for scaffolds > " + MIN_SCAFFOLD + " from proportion of coverage = " + n50scaffold +
                     ";\nTotal Scaffold Length from proportion of coverage = " + myFormatter.format(pSC.getTotalAllScaffLength()) +
                     ";\nN50 Scaffold Length for all scaffolds from proportion of coverage = " + n50AllScaffold;
            }
            if (!CONTIG_GAPS_FILE.getName().contentEquals("temp.tmp")) {
               printHistogram(countContigsInScaffolds, COUNT_CONTIGS_IN_SCAFFOLDS_FILE);
            }
         }


         if (!CONTIGS_FILE.getName().contentEquals("temp.tmp")) {
            printCoordinates(CONTIGS_FILE, FileType.bed, contigs);
         }
      }

      if (!CONTIG_GAPS_FILE.getName().contentEquals("temp.tmp")) {
         printHistogram(contigGaps, CONTIG_GAPS_FILE);
      }

      if (!CONTIG_LENGTHS_FILE.getName().contentEquals("temp.tmp")) {
         printHistogram(contigLengths, CONTIG_LENGTHS_FILE);
      }

      if (!COVERAGE_FILE.getName().contentEquals("temp.tmp")) {
         printHistogram(coverage, COVERAGE_FILE);
      }


      try {
         //Start writing to the output stream
         writer.write(outString);
         writer.newLine();

      }
      catch (FileNotFoundException ex) {
         ex.printStackTrace();
      }
      catch (IOException ex) {
         ex.printStackTrace();
      }
      finally {
         closeStreams(writer);

      }

   }

   private void closeStreams(BufferedWriter stream) {
      //Close the BufferedWriter
      try {
         if (stream != null) {
            stream.flush();
            stream.close();
         }

      }
      catch (IOException ex) {
         ex.printStackTrace();
      }

   }

   private String getN50(HashMap<Integer, Long> hmap) {
      Integer test = 50;
      Long n50 = 0L;
      while (hmap.get(test) == null && test < 100) {
         test++;
      }
      if (test == 50) {
         n50 = hmap.get(test);
      }
      else {
         int lower = 50;
         while (hmap.get(lower) == null && lower > 0) {
            lower--;
         }
         n50 = hmap.get(lower) + ((50 - lower) / (test - lower)) * (hmap.get(test) - hmap.get(lower));
      }
      if (n50 != 0) {
         return "50th Percentile = " + myFormatter.format(n50);
      }
      else {
         return " = Not found; Graph length " + hmap.size();
      }
   }

   private void printPercentile(int percent, long start, long end, int intervalCount) {
      String instance = myFormatter.format(start);
      String cumulative = myFormatter.format(end);
      String interCount = myFormatter.format(intervalCount);
      try {
         writer.write("Percentile\t" + inputFilename + "\t" + percent + "\t" + instance + "\t" + cumulative + "\t" + interCount);
         writer.newLine();
      }
      catch (IOException ex) {
         ex.printStackTrace();
      }
   }

   private void printArray(ArrayList<Integer> array, String arrayName, Integer minValue) {
      for (Integer value : array) {
         if (value > minValue) {
            String instance = myFormatter.format(value);
            try {
               writer.write("Percentile\t" + inputFilename + "\t" + arrayName + "\t" + instance);
               writer.newLine();
            }
            catch (IOException ex) {
               ex.printStackTrace();
            }

         }

      }
   }

   private void printHeader(String header) {
      try {
         writer.write(header);
         writer.newLine();
      }
      catch (IOException ex) {
         ex.printStackTrace();
      }

   }

   private void printHistogram(HashMap<Long, Integer> lengths, File file) {

      IoUtil.assertFileIsWritable(file);
      BufferedWriter hist_writer = IoUtil.openFileForBufferedWriting(file, false);


      ArrayList<Long> keys = new ArrayList<Long>();
      keys.addAll(lengths.keySet());
      // and sorting it
      Collections.sort(keys);

      for (Long i : keys) {
         try {
            hist_writer.write(i.toString() + "\t" + lengths.get(i).toString());
            hist_writer.newLine();
         }
         catch (IOException ex) {
            ex.printStackTrace();
         }

      }
      closeStreams(hist_writer);
   }

   private void printReadPairs() {
      String filename = inputFilename + "_ReadPairCoordinates.txt";
      filename = filename.replaceAll("\\s", "");
      File readPairFile = new File(filename);
      printCoordinates(readPairFile, FileType.bed, readPairs);
   }

   //Print arraysList of Coordinates can be scaffolds or contigs
   private void printCoordinates(File file, FileType type, ArrayList<Coordinates> coords) {
      IoUtil.assertFileIsWritable(file);
      BufferedWriter bWriter = IoUtil.openFileForBufferedWriting(file, false);
      for (Coordinates coord : coords) {
         try {
            if (type.equals(FileType.bed)) {
               bWriter.write(coord.getIntervalBed());
            }
            if (type.equals(FileType.gff3)) {
               bWriter.write(coord.getIntervalGff3());
            }
            if (type.equals(FileType.interval_list)) {
               bWriter.write(coord.getInterval());
            }
            bWriter.newLine();
         }
         catch (IOException ex) {
            ex.printStackTrace();
         }
      }
      closeStreams(bWriter);

   }

   private class Coordinates implements Comparable<Coordinates> {

      private long start;
      private long end;
      private String chr;
      private String recName;
      private Double coverage = -1.0;
      private final Double index;//sum of start position and a fractional offset ot make multiple different reads startign at same position unique
      private SofaType type;

      public Coordinates(String recName, String chr, long start, long end, Double index) {
         this.index = index;
         this.recName = recName;
         this.start = start;
         this.end = end;
         this.chr = chr;

      }

      public Coordinates(String recName, String chr, long start, long end) {
         this.recName = recName;
         this.start = start;
         this.end = end;
         this.chr = chr;
         this.index = start + 0.01;
      }

      public Coordinates(String recordName, String chr, int start) {
         this.start = start;
         this.chr = chr;
         this.index = start + 0.01;
      }

      public void setEnd(int end) {
         this.end = end;
      }

      public String getChr() {
         return chr;
      }

      public Double getIndex() {
         return index;
      }

      public void setCoverage(Double cov) {
         this.coverage = cov;
      }

      public Double getCoverage() {
         return coverage;
      }

      public void setType(SofaType type) {
         this.type = type;
      }

      public SofaType getType() {
         return type;
      }

      public long getEnd() {
         return end;
      }

      public long getStart() {
         return start;
      }

      public String getInterval() {
         return recName + ":" + chr + ":" + start + "-" + end;
      }

      public String getIntervalBed() {
         return chr + "\t" + start + "\t" + end;
      }

      public String getIntervalGff3() {
         String source = "PicardTools ContigMetrics.jar";
         return chr + "\t" + source + "\t" + type + "\t" + start + "\t" + end + "\t.\t.\t.\tcoverage=" + floatFormat.format(coverage);
      }

      public int compareTo(Coordinates co) {
         int cmp = index.compareTo(co.getIndex());
         return cmp;
      }
   }

   //enum class for type of GFF3
   // se http://www.sequenceontology.org/gff3.shtml
   //http://www.sequenceontology.org/browser/current_cvs/term/SO:0000148 supercontig (scaffold)
   //http://www.sequenceontology.org/browser/current_cvs/term/SO:0000149 contig
   public enum SofaType {

      contig, scaffold;
   };

   public enum FileType {

      gff3, bed, interval_list;
   };

//helper class for scaffold metrics
   private class ScaffoldMetrics {

      private long contCount = 0;
      private int scaffCount = 0;
      private long totalScaffLength = 0;
      private int allScaffCount = 0;
      private long totalAllScaffLength = 0;
      private ArrayList<Long> allScaffLengths;
      private ArrayList<Long> scaffLengths;
      private Double stdev;
      private Double maxGap;
      private ArrayList<Coordinates> scaffolds;

      public ScaffoldMetrics(ArrayList<Coordinates> scaffolds, Double maxGap) {
         this.scaffolds = scaffolds;
         this.maxGap = maxGap;

         initialise();
      }

      protected void initialise() {
         scaffCount = 0;
         allScaffCount = 0;
         totalScaffLength = 0;
         scaffLengths = new ArrayList<Long>();
         allScaffLengths = new ArrayList<Long>();
         for (Coordinates coords : scaffolds) {
            long scaffLength = coords.end - coords.start;
            if (scaffLength > MIN_SCAFFOLD) {
               scaffCount++;
               totalScaffLength += scaffLength;
               totalAllScaffLength += scaffLength;
               scaffLengths.add(scaffLength);
               allScaffLengths.add(scaffLength);
               allScaffCount++;
            }
            else {
               totalAllScaffLength += scaffLength;
               allScaffLengths.add(scaffLength);
               allScaffCount++;
            }
         }


      }

      public Double getMaxGap() {
         return maxGap;
      }

      public void setStdev(Double stdev) {
         this.stdev = stdev;
      }

      public ArrayList<Coordinates> getScaffolds() {
         return scaffolds;
      }

      public void setScaffolds(ArrayList<Coordinates> scaffolds) {
         this.scaffolds = scaffolds;
         initialise();
      }

      public int getAllScaffCount() {
         return allScaffCount;
      }

      public ArrayList<Long> getAllScaffLengths() {
         return allScaffLengths;
      }

      public long getTotalAllScaffLength() {
         return totalAllScaffLength;
      }

      public void setContigCount(long contCount) {
         this.contCount = contCount;
      }

      public long getContCount() {
         return contCount;
      }

      public Double getStdev() {
         return stdev;
      }

      public int getScaffCount() {
         return scaffCount;
      }

      public ArrayList<Long> getScaffLengths() {

         return scaffLengths;
      }

      public long getTotalScaffLength() {
         return totalScaffLength;
      }
   }
}
