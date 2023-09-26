# perlotSeeds
Periodic lossless ternary seeds of maximum weight

<nav>
  <ul>
    <li><a href="#link_intro">Introduction</a></li>
    <li><a href="#link_binaryLevel">periodicBinaryBlockLevel: binary blocks of less than maximum weight</a></li>
    <li><a href="#link_ternaryBlock">periodicTernaryBlocks: ternary blocks of maximum weight</a></li>	  
 <li><a href="#link_bin2text">convertBin2Text: convert ternary blocks in binary format into a text</a></li>
    <li><a href="#link_maxWeight">ternarySeedMaxWeight: ternary seeds of maximum weight</a></li>
  </ul>
  </nav>

<h2 id="link_intro">Introduction</h2>
We consider sequences of symbols <tt>A</tt>, <tt>C</tt>, <tt>G</tt> and <tt>T</tt>. In practical applications we have a long <i>reference</i> sequence (upto billions of symbols) and a short sequence (<i>read</i>, often 50-300 symbols) obtained experiemntally. The read is a small chunk of an unknown sequeunce, which is in some sense similar to the reference sequence. Therefore we expect these two long sequences to have similar subsequnces. There will be also some deviations, like SNPs (single-nucleotide polymorphisms) when some symbols are replaced by other symbols or insertions/deliations when some symbols are added or removed. 

The simplest approach to form the unknown long sequence is to pre-align its chunks (reads) with respect to the known refernce sequence, then perform full scale comparision to take into account possible SNPs and insersions/deletions. For the pre-aligning step we usually create a library of pairs (position within the reference sequence, the sequence of symbols) and sort all records by "sequence of symbols" value. So, if a read contains one of the subseqeunces from the library's records, then we can find the corresponding positions within the reference sequnce and thus pre-align the read. Of course, this approach only works if all elements of the chunks are the same.

For example, we have two sequences <tt>TTGGAGATCG</tt> and <tt>TAGGTGCTCG</tt> (of length 10). We compare the corresponding elements of the sequences and get 1 if there is a match and 0 if there is a mismatch.

<table>
  <tr><th>Sequence 1</th><th><tt>TTGGAGATCG</tt></th></tr>
  <tr><th>Sequence 2</th><th><tt>TAGGTGCTCG</tt></th></tr>
  <tr><th>Match</th><th><tt>1011010111</tt></th></tr>
</table>

We see that if we know the poistion of the first sequence <tt>TTGGAGATCG</tt> (for example, position 16 with respect to the reference <tt>ACGACAACCTTGTCGTTGGAGATCGGAAGAGCACACGTCTGAAC</tt>), then we may pre-align another string containing <tt>TTGGAGATCG</tt>. However, it is not possible to pre-align a string containing <tt>TAGGTGCTCG</tt>, since the library of records does not contain this chunk.

Reads often SNPs, so it is cruicial to deal with possible mismatches. One of the approaches is to use <i>seeds</i>, i.e. a sequence of 0 and 1 elements. Let there be two seqeunce of symbols and a seed of the same length. When an element of the seed is 1, then we compare corresponding elements of two symbol sequences, otherwise we ignore possible deviations. 

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
So, if we form a library of record based on <tt>1010010101</tt>, then one record is (16, <tt>TGGTG</tt>) and a read contining the second sequence <tt>TAGGTGCTCG</tt> can be pre-aligned with respect to <tt>ACGACAACCTTGTCGTTGGAGATCGGAAGAGCACACGTCTGAAC</tt>.

Ideally, we should find seeds of large weight, since by increasing the weight by one we reducing the number of candidate positions to be chekced in 4 times (aassuming that the chance to have any of the symbols <tt>A</tt>, <tt>C</tt>, <tt>G</tt> and <tt>T</tt> is the same). Codes to generate such seed can be found in <a href="https://github.com/vtman/PerFSeeB">https://github.com/vtman/PerFSeeB</a>.

Standard seeds are <i>binary</i> when there are only two states (<b><tt>1</tt></b> or <b><tt>#</tt></b> for <i>match</i> and <b><tt>0</tt></b> or <b><tt>&#95;</tt></b> for ``don't care symbol''). In genetics, the chance to have a <b>transition</b> mutation (<tt>A</tt> &harr; <tt>G</tt> or <tt>C</tt> &harr; <tt>T</tt>) is often twice higher than a <b>transversion</b> mutation (<tt>A</tt> &harr; <tt>C</tt>, <tt>A</tt> &harr; <tt>T</tt>, <tt>G</tt> &harr; <tt>C</tt>, <tt>G</tt> &harr; <tt>T</tt>). Transition-constrained seeds use ternary alphabet {<b><tt>#</tt></b>, <b><tt>@</tt></b>, <b><tt>&#95;</tt></b>} where <b><tt>@</tt></b> is for a match or a transition mismatch.

We construct ternary seeds as periodic seeds when there is a whole number of repeatitions of a string followed by a remainder. 


<h2 id="link_binaryLevel">periodicBinaryBlockLevel: binary blocks of less than maximum weight</h2>

In PerFSeeB project we have gnerated periodic binary blocks of maximum weight for a given number of mismatches. We may reuse results obtained in that project and also generate binary blocks that have less than the maximum weight.

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

In most case it is enough to use level = 0 (85%), level = 2 only requires in a couple cases, the rest is for level = 1.

Output files can be downloaded from <a href="zenodo.com">Zenodo</a> and example output fiels are <a href="https://github.com/vtman/perlotSeeds/tree/main/ExampleBinaryBlocks">ExampleBinaryBlocks</a>.

Output files are in binary format. For a given block size <b>B</b> we find the smallest number <b>N</b> such that <b>B &#8804;8N</b>. So, for block of length 30, we need 4 bytes to store. A hundred blocks requires 100*4= 400 bytes.

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

Output files for ternary blocks are also in binary format. Each element of a block requires 2 bits (<tt>_</tt> = <tt>0<tt> = <tt>00<tt>, <tt>#</tt> = = <tt>1<tt> = <tt>01<tt>, <tt>@</tt> = = <tt>2<tt> = <tt>10<tt>). So, a block of length 30 requires 8 bytes.

<h2 id="link_bin2text">convertBin2Text: convert ternary blocks in binary format into a text</h2>

<h3>Parameters</h3>

<ol>
  <li>Input folder (ternary blocks in binary format)</li>
  <li>Output folder</li>
  <li>Block size</li>
  <li>Number of mismatches (transition)</li>
  <li>Number of mismatches (transversion)</li>
<li>Level</li>
</ol>

<tt>convertBin2Text.exe E:\Temp2\perlotSeeds\ternaryBlocks E:\Temp2\perlotSeeds\ternaryBlocksText 30 2 3</tt>


<h2 id="link_maxWeight">ternarySeedMaxWeight: ternary seeds of maximum weight</h2>

For practical applications, it is better to use seeds of maximum weight. We have four letters <tt>A</tt>, <tt>C</tt>, <tt>G</tt> and <tt>T</tt>. Let us assume that their chance to be in a sequence is the same and does not depend on neighbouring symbols (these assumptions are not completely true). So, if a chance to find a pattern of length k within a reference sequence is P<sub>k</sub>, then the chance to find a pattern of length (k+1) is P<sub>k</sub>/4. Therefore seeds of higher weight allow us to process 4 times fewer candidate positions within a reference sequence.

For a given length of reads and a number of mismatches, there may be several seeds of maximum weight. For example, for r=45, m=8 and w=6 we get several valid seeds including <tt>1111011</tt>, <tt>1011110001</tt>, <tt>11011000101</tt>, <tt>1001001000001000000001001</tt>, <tt>1100000001000000001100000001</tt>. This means that if we choose one of the shortest seeds (<tt>1111011</tt>) of length 7, then we need to consider 45-7+1=39 chunks of a read (of length 45) and use them to find corresponding candidate positions within a reference seqeunce. However, if we use seed <tt>1100000001000000001100000001</tt> of length 28, the number of chunks becomes 45-28+1=18 (almost 2 times less). So, using longest seeds among seeds of maximum weight can reduce processing times.

Of course, there may several lengths of reads when w is the maximum weight. If a seed is valid for a read of length r, then it is also valid for a read of length (r+1). Therefore, while there may be several seeds valid for various lengths of reads, we pick up only those valid for the shortest lengths. For around 80% of seeds obtained using the iterative procedure we may see that the best seeds (longest seeds of maximum weight valid for shortest reads) have a periodic structure: an integer number n<sub>b</sub> of blocks of length <i>T</i> and a "remainder" (first n<sub>d</sub> elements of the block), so the total length is n<sub>s</sub> = n<sub>b</sub> T + n<sub>d</sub>, and the following formula is valid

r = n<sub>s</sub> + T - 1

For example, seed <tt>1110100000000111010000000011101</tt> is valid for m=4, r=43 and can be split up as
<table>
	<tr><th><tt>1110100000000</tt></th><th><tt>1110100000000</tt></th><th><tt>11101</tt></th></tr>
</table>
We seed that T=13, n<sub>b</sub>=2, n<sub>d</sub>=5, so n<sub>s</sub>=31 and 43 = 31 + 13 - 1.

We try to find possible blocks such that we can form seeds of the given structure.




