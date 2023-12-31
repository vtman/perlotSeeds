# perlotSeeds
Periodic lossless ternary seeds of maximum weight

<nav>
  <ul>
    <li><a href="#link_intro">Introduction</a></li>
    <li><a href="#link_binaryLevel">periodicBinaryBlockLevel: binary blocks of less than maximum weight</a></li>
    <li><a href="#link_ternaryBlock">periodicTernaryBlocks: ternary blocks of maximum weight</a></li>	  
 <li><a href="#link_bin2text">convertBin2Text: convert ternary blocks in binary format into a text</a></li>
    <li><a href="#link_maxWeight">ternarySeedMaxWeight: ternary seeds of maximum weight</a></li>
    <li><a href="#link_data">Linked data</a></li>
  </ul>
  </nav>

<h2 id="link_intro">Introduction</h2>
We consider sequences of symbols <tt>A</tt>, <tt>C</tt>, <tt>G</tt> and <tt>T</tt>. In practical applications, we have a long <i>reference</i> sequence (up to billions of symbols) and a short sequence (<i>read</i>, often 50-300 symbols) obtained experimentally. The read is a small chunk of an unknown sequence, similar to the reference sequence. Therefore, we expect these two long sequences to have similar subsequences. There will also be some deviations, like SNPs (single-nucleotide polymorphisms), when some symbols are replaced by other symbols or insertions/deletions when some symbols are added or removed. 

The most straightforward approach to form the unknown long sequence is to pre-align its chunks (reads) to the known reference sequence and then perform a full-scale comparison to consider possible SNPs and insertions/deletions. For the pre-aligning step, we usually create a library of pairs (position within the reference sequence, the sequence of symbols) and sort all records by "sequence of symbols" value. So, if a read contains one of the subsequences from the library's records, we can find the corresponding positions within the reference sequence and thus pre-align the read. Of course, this approach only works if all elements of the chunks are the same.

For example, we have two sequences <tt>TTGGAGATCG</tt> and <tt>TAGGTGCTCG</tt> (of length 10). We compare the corresponding elements of the sequences and get 1 if there is a match and 0 if there is a mismatch.

<table>
  <tr><th>Sequence 1</th><th><tt>TTGGAGATCG</tt></th></tr>
  <tr><th>Sequence 2</th><th><tt>TAGGTGCTCG</tt></th></tr>
  <tr><th>Match</th><th><tt>1011010111</tt></th></tr>
</table>

We see that if we know the position of the first sequence <tt>TTGGAGATCG</tt> (for example, position 16 to the reference <tt>ACGACAACCTTGTCGTTGGAGATCGGAAGAGCACACGTCTGAAC</tt>), then we may pre-align another string containing <tt>TTGGAGATCG</tt>. However, it is not possible to pre-align a string containing <tt>TAGGTGCTCG</tt> since the library of records does not contain this chunk.

Reads often SNPs, so it is crucial to deal with possible mismatches. One of the approaches is to use <i>seeds</i>, i.e. a sequence of 0 and 1 elements. Let two sequences of symbols and a seed of the same length. When an element of the seed is 1, we compare the corresponding elements of two symbol sequences; otherwise, we ignore possible deviations. 

If we use seed <tt>1010101010</tt>, then we only need to compare symbols at odd positions, five in total. Three symbols match in the sequences. 
<table>
  <tr><th>Sequence 1</th><th><tt>TTGGAGATCG</tt></th></tr>
  <tr><th>Sequence 2</th><th><tt>TAGGTGCTCG</tt></th></tr>
  <tr><th>Seed</th>      <th><tt>1010101010</tt></th></tr>
  <tr><th>Match</th>     <th><tt>1_1_0_0_1_</tt></th></tr>
</table>

If we use seed <tt>1010010101</tt> (length is 10, weight is 5, weight is the number of 1-elements), then we have all five symbols match.
<table>
  <tr><th>Sequence 1</th><th><tt>TTGGAGATCG</tt></th></tr>
  <tr><th>Sequence 2</th><th><tt>TAGGTGCTCG</tt></th></tr>
  <tr><th>Seed</th>      <th><tt>1010010101</tt></th></tr>
  <tr><th>Match</th>     <th><tt>1_1__1_1_1</tt></th></tr>
</table>
So, if we form a library of records based on <tt>1010010101</tt>, then one record is (16, <tt>TGGTG</tt>) and a read containing the second sequence <tt>TAGGTGCTCG</tt> can be pre-aligned to <tt>ACGACAACCTTGTCGTTGGAGATCGGAAGAGCACACGTCTGAAC</tt>.

Ideally, we should find seeds of considerable weight since by increasing the weight by one, we reduce the number of candidate positions to be checked in 4 times (assuming that the chance to have any of the symbols <tt>A</tt>, <tt>C</tt>, <tt>G</tt> and <tt>T</tt> is the same). Codes to generate such seed can be found in <a href="https://github.com/vtman/PerFSeeB">https://github.com/vtman/PerFSeeB</a>.

Standard seeds are <i>binary</i> when there are only two states (<b><tt>1</tt></b> or <b><tt>#</tt></b> for <i>match</i> and <b><tt>0</tt></b> or <b><tt>&#95;</tt></b> for ``don't care symbol''). In genetics, the chance to have a <b>transition</b> mutation (<tt>A</tt> &harr; <tt>G</tt> or <tt>C</tt> &harr; <tt>T</tt>) is often twice higher than a <b>transversion</b> mutation (<tt>A</tt> &harr; <tt>C</tt>, <tt>A</tt> &harr; <tt>T</tt>, <tt>G</tt> &harr; <tt>C</tt>, <tt>G</tt> &harr; <tt>T</tt>). Transition-constrained seeds use ternary alphabet {<b><tt>#</tt></b>, <b><tt>@</tt></b>, <b><tt>&#95;</tt></b>} where <b><tt>@</tt></b> is for a match or a transition mismatch.

We construct ternary seeds as periodic seeds when there is a whole number of repetitions of a string followed by a remainder. 


<h2 id="link_binaryLevel">periodicBinaryBlockLevel: binary blocks of less than maximum weight</h2>

In the PerFSeeB project, we have generated periodic binary blocks of the maximum weight for a given number of mismatches. We may reuse results obtained in that project and generate binary blocks with less than the maximum weight.

<h3>Parameters</h3>

<ol>
  <li>Input folder (binary blocks obtained in PerFSeeB project, can be found in <a href="https://github.com/vtman/perlotSeeds/tree/main/binaryBlocksText">binaryBlocksText</a>)</li>
  <li>Output folder</li>
  <li>Number of mismatches (from 2 to 9)</li>
  <li>Size of blocks (from 3 to 50)</li>
  <li>Minimum level (0 when maximum weight)</li>
  <li>Maximum level</li>
</ol>

<tt>periodicBinaryBlockLevel.exe E:\Temp2\perlotSeeds\binaryBlocksText E:\Temp2\perlotSeeds\binaryBlocks 4 24 0 2</tt>

In most cases, it is enough to use level = 0 (85%); level = 2 is only required in a couple of cases, and the rest is for level = 1.

Output files can be downloaded from <a href="https://zenodo.org/record/8370909">Zenodo</a>, and example output files are <a href="https://github.com/vtman/perlotSeeds/tree/main/ExampleBinaryBlocks">ExampleBinaryBlocks</a>.

Output files are in binary format. For a given block size <b>B</b>, we find the smallest number <b>N</b> such that <b>B &#8804;8N</b>. So, for a block of length 30, we need 4 bytes to store. A hundred blocks require 100*4= 400 bytes.

<h2 id="link_ternaryBlock">periodicTernaryBlocks: ternary blocks of maximum weight</h2>


<h3>Parameters</h3>

<ol>
  <li>Input folder (files for binary blocks, like in <a href="https://github.com/vtman/perlotSeeds/tree/main/ExampleBinaryBlocks">ExampleBinaryBlocks</a>)</li>
  <li>Output folder</li>
  <li>Block size</li>
  <li>Number of mismatches (transition)</li>
  <li>Number of mismatches (transversion)</li>
<li>Level</li>
</ol>

<tt>periodicTernaryBlocks.exe E:\Temp2\perlotSeeds\binaryBlocks E:\Temp2\perlotSeeds\ternaryBlocks 30 2 3 0</tt>

Output files for ternary blocks are also in binary format. Each element of a block requires 2 bits (<tt>&#95;</tt> = <tt>0</tt> = <tt>00</tt>, <tt>#</tt> = <tt>1</tt> = <tt>01</tt>, <tt>@</tt> = <tt>2</tt> = <tt>10</tt>). So, a block of length 30 requires 8 bytes.

Output files can be downloaded from <a href="https://zenodo.org/record/8370909">Zenodo</a>, and example output files are <a href="https://github.com/vtman/perlotSeeds/tree/main/ExampleTernaryBlocks">ExampleTernaryBlocks</a>. Some files are huge (tens of gigabytes). 

<h2 id="link_bin2text">convertBin2Text: convert ternary blocks in binary format into a text</h2>

<h3>Parameters</h3>

<ol>
  <li>Input folder (ternary blocks in binary format)</li>
  <li>Output folder</li>
  <li>Block size</li>
  <li>Number of mismatches (transition)</li>
  <li>Number of mismatches (transversion)</li>
</ol>

<tt>convertBin2Text.exe E:\Temp2\perlotSeeds\ternaryBlocks E:\Temp2\perlotSeeds\ternaryBlocksText 30 2 3</tt>


<h2 id="link_maxWeight">ternarySeedMaxWeight: ternary seeds of maximum weight</h2>

<h3>Parameters</h3>

<ol>
  <li>Input folder (ternary blocks in binary format)</li>
  <li>Output folder</li>
  <li>Number of mismatches (transition)</li>
  <li>Number of mismatches (transversion)</li>
  <li>Block size (minimum)</li>
  <li>Block size (maximum)</li>
  <li>Read length (minimum)</li>
  <li>Read length (maximum)</li>
</ol>

<tt>ternarySeedMaxWeight.exe E:\Temp2\TSeeds\Zen\TernaryBlocks\T3V4 E:\Temp2\perlotSeeds\ternarySeeds 3 4 3 50 50 80</tt>

<h2 id="link_data">Linked data</h2>

Binary, ternary blocks and best ternary seeds for reads of length 30 to 500 can be found in Zenodo archives (<a href="https://zenodo.org/record/8395215">10.5281/zenodo.8395215</a> and <a href="https://zenodo.org/record/8395813">10.5281/zenodo.8395813</a>).
