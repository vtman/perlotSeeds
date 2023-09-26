# perlotSeeds
Periodic lossless ternary seeds of maximum weight

<nav>
  <ul>
    <li><a href="#link_intro">Introduction</a></li>
    <li><a href="#link_binaryLevel">periodicBinaryBlockLevel: binary blocks of less than maximum weight</a></li>
    <li><a href="#link_ternaryBlock">periodicTernaryBlocks: ternary blocks of maximum weight</a></li>	  
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

Output files are in binary format. For a given block size <b>B<b> we find the smallest number <b>N</b> such that <b>B &#8804;8N</b>. So, for block of length 30, we need 4 bytes to store. A hundred blocks requires 100*4= 400 bytes.

<h2 id="link_ternaryBlock">periodicTernaryBlocks: ternary blocks of maximum weight</h2>


<h3>Parameters</h3>

<ol>
  <li>Input folder (files for binary blocks, like in <a href="https://github.com/vtman/perlotSeeds/tree/main/ExampleBinaryBlocks">ExampleBinaryBlocks</a>)</li>
  <li>Outout folder</li>
  <li>Block size</li>
  <li>Number of mismatches (transition)</li>
  <li>Number of mismatches (transversion)</li>
<li>Level</li>
</ol>


Since there are usually a lot of spaced seeds generated in this way, we try to report only seeds of large weights. Therefore we specify the minimum weight required for a seed to be reported.

<tt>iterSeed.exe C:\MyFolder 15 2 6</tt>

As an output we get the following spaced seeds: <tt>111010011</tt>, <tt>111001011</tt>, <tt>110100111</tt>, <tt>110010111</tt>, <tt>1101001101</tt>, <tt>1011001011</tt>. These seeds have the maximum weight possible (6). Note that if there is a valid seed, then its flipped version should also be in the list, i.e. seeds <tt>110100111</tt>, <tt>110010111</tt> are both in the list. There may also be longer seeds of smaller weight. For example, the maximum length for seeds of weight 6 is 10 (seeds <tt>1101001101</tt> and <tt>1011001011</tt>), however for seeds of weight 5 we get <tt>1001001001001</tt> of length 13.

We have generated spaced seeds for a number of mismatches from 2 to 8 and reads' lengths from 10 to 50 (different ranges for different numbers of mismatches). They are in <b>iterSeed</b> folder.

<h2 id="link_maxWeight">Seeds of maximum weight</h2>

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

<h2 id="link_periodicBlock">periodicBlock: Periodic blocks</h2>

If the formula above is valid for a periodic seed, then to validate the seed is enough to validate its periodic block. For example, we have seed <tt>110010111001011100101110010111</tt>. We create 7 rows such that

<table>
	<tr><th><tt>110010111001011100101110010111000000</tt></th></tr>
	<tr><th><tt>011001011100101110010111001011100000</tt></th></tr>
	<tr><th><tt>001100101110010111001011100101110000</tt></th></tr>
	<tr><th><tt>000110010111001011100101110010111000</tt></th></tr>
	<tr><th><tt>000011001011100101110010111001011100</tt></th></tr>
	<tr><th><tt>000001100101110010111001011100101110</tt></th></tr>
	<tr><th><tt>000000110010111001011100101110010111</tt></th></tr>
</table>

We split the rows to form
<table>
	<tr><th><tt>1100101</tt></th><th><tt>1100101</tt></th><th><tt>1100101</tt></th><th><tt>1100101</tt></th><th><tt>11000000</tt></th></tr>
	<tr><th><tt>0110010</tt></th><th><tt>1110010</tt></th><th><tt>1110010</tt></th><th><tt>1110010</tt></th><th><tt>11100000</tt></th></tr>
	<tr><th><tt>0011001</tt></th><th><tt>0111001</tt></th><th><tt>0111001</tt></th><th><tt>0111001</tt></th><th><tt>01110000</tt></th></tr>
	<tr><th><tt>0001100</tt></th><th><tt>1011100</tt></th><th><tt>1011100</tt></th><th><tt>1011100</tt></th><th><tt>10111000</tt></th></tr>
	<tr><th><tt>0000110</tt></th><th><tt>0101110</tt></th><th><tt>0101110</tt></th><th><tt>0101110</tt></th><th><tt>01011100</tt></th></tr>
	<tr><th><tt>0000011</tt></th><th><tt>0010111</tt></th><th><tt>0010111</tt></th><th><tt>0010111</tt></th><th><tt>00101110</tt></th></tr>
	<tr><th><tt>0000001</tt></th><th><tt>1001011</tt></th><th><tt>1001011</tt></th><th><tt>1001011</tt></th><th><tt>10010111</tt></th></tr>
</table>

To validate the seed we need to choose any m arbitrary columns. For those columns, there should be a row where all these columns have 0-elements. Therefore if choose columns within the second "big" column (of 7 standard columns), then the validation procedure should be true for it. All other columns outside this "big"  column either have identical columns within the "big" one or have more 0-elements than corresponding columns within the "big" one. Therefore we only need to check the "big" column, which is made by the cyclic shift operation.
<table>
	<tr><th><tt>1100101</tt></th></tr>
	<tr><th><tt>1110010</tt></th></tr>
	<tr><th><tt>0111001</tt></th></tr>
	<tr><th><tt>1011100</tt></th></tr>
	<tr><th><tt>0101110</tt></th></tr>
	<tr><th><tt>0010111</tt></th></tr>
	<tr><th><tt>1001011</tt></th></tr>
</table>

<h3>Parameters</h3>

<ol>
  <li>Output folder</li>
  <li>Number of mismatches</li>
  <li>Length of a block</li>
  <li>Initial number of 0-elements</li>
  <li>Size of chunks for parallelisation</li>
  <li>Number of threads</li>
  <li>(optional) omit first N candidates (default = -1)</li>
  <li>(optional) Initial string </li>
</ol>

<tt>periodicBlock.exe D:\PerFSeeB\outTest 3 30 2 100000 6</tt>

Our goal is to find blocks of maximum weight (maximum number of 1-elements). Therefore we start we a small number of 0-elements and try to find any valid blocks. If there are no such blocks, then we increase the number of 0-elements in a block by one and repeat the procedure. For short blocks (< 35) the procedure is fast, so the number can be any small (and must be positive). However, for big blocks, it may be reasonable to start with some number, e.g. using numbers found for shorter blocks. However, one needs to take into account that while the number of 0-elements tends to increase with the block size, some blocks may have fewer zeros than shorter blocks.

The code pre-generates a list of candidate blocks, then all those candidates are processed (validated) in parallel. Therefore a user must specify the number of candidates to be pre-generated and the number of threads to be used. The number of candidates is good to set to 1000000, the number of threads is usually the number of cores in a CPU (in any case the code will check the number of threads available).

Processing for large blocks may take hours/days, so in some cases when the code is stopped one may use the last found block as the initial block. Or one may specify the number of candidate blocks to be skipped (generation of candidate blocks is relatively fast compared to the time required for validation). However, these two parameters are optional and may not be used in a normal situation.

We aim for blocks of maximum weight, therefore once at least one block of a given weight is found, we do not generate blocks of smaller weights. In principle, if a user needs other blocks, then it is possible to modify a code to specify the exit procedure or use another number of 0-elements. Note that the number of blocks can be large in the case of non-maximum weights.

Blocks generated by the above code are saved in <b>perBlock</b> folder.

<h2 id="link_bestPerSeed">bestPerSeed: Finding best periodic seeds</h2>

Blocks found above can be used to generate spaced seeds of high weight. For a given lengths of reads we try various periodic blocks, so that r = n<sub>S</sub> + T - 1. For the first step, we identify the maximum weight of spaced seeds. For the second step, we write down those seeds into files. 

<h3>Parameters</h3>

<ol>
  <li>Input folder</li>
  <li>Output folder</li>
  <li>Number of mismatches</li>
  <li>File index (first)</li>
  <li>File index (last)</li>
  <li>Read length (min)</li>
  <li>Read length (max)</li>
</ol>

<tt>bestPerSeed.exe D:\PerFSeeB\perBlocks D:\PerFSeeB\outTest\bestS 5 6 50 20 200</tt>

Input files are those generated by <i>periodicBlock</i>. The code uses the same template of file names. File indices are the minimum/maximum lengths of blocks to be used. Output files are generated for given lengths of reads. Each output file contains the maximum weight, size of a periodic block T, the number of these periodic blocks n<sub>b</sub>, "remainder" n<sub>d</sub>, the periodic block and the final seed. There are usually several seeds in each file. For short reads, there may be no files generated (empty files).

There is also <b>res_#.txt</b> in the output folder containing information about read's length, seed weight, number of sizes for periodic blocks that give us the highest weight and (after | separator) the list of those sizes. Note that the code generates seeds based on the main formula. Therefore in some cases, seeds of higher weight may be found for shorter reads. Thus it is important to check <b>res_#.txt</b> files to identify those lengths.

<h2  id="link_bestLaTeX">bestSeedsLaTeX</h2>
For a given number of mismatches and seed's weight, we find the minimum length of reads that can be used and report those best seeds. The code has the same parameters and generates an itemised list for a LaTeX file. By default, we restrict the weights (between 16 and 320), it can be changed in the code (countLoc variable).

PDF files for the latest seeds are in <b>bestSeeds</b> folder.


<h2 id="link_tools">Tools</h2>

To deal with sequences we convert them into a binary file format. We use two files formats:
<ol>
  <li><b>Standard</b>. We assume that we have only four symbols (<tt>A</tt>, <tt>C</tt>, <tt>G</tt>, <tt>T</tt>). Then each each letter can be coded as a 2-bit symbol (<tt>A=0</tt>, <tt>C=1</tt>, <tt>G=2</tt>, <tt>T=3</tt>). So each for each sequence we store its length <tt>L</tt> (usually in a separete file), the find the smallest number <tt>K</tt> such that <tt>L &leq; 4*K</tt> and generate a binary sequence of <tt>8*K</tt> bits (or, equivalently, <tt>K</tt> bytes). We pad the original sequence with <tt>A</tt> symbols if it is needed. For some problems we require <tt>L &leq; 32*K</tt>.</li>
  
  <li><b>m128</b>. If <tt>N</tt>symbols are to be taken into account or we need to count the number of same symbols for two seqeunces, then we use a different file format. We split the sequence into groups of 32 symbols. We write this 32-symbol sequence as a 128-bit number. An <i>i</i>-th bit from the first 32 bits is 1 if the <i>i</i>-th symbol in the sequence is <tt>A</tt>, otherwise it is 0. The next 32 bits are for symbol <tt>C</tt>, then for symbol <tt>G</tt> and <tt>T</tt> respectively. If there is a symbol <tt>N</tt>, then all bits in <tt>A</tt>, <tt>C</tt>, <tt>G</tt> and <tt>T</tt> arrays are 0.</li>
</ol>

<hr>
<b>Example 1</b>

Let us have the following sequence of 32 symbols: <tt>CATAGNCACGTGATCCTAGNCATGTTACCTGT</tt>.

<table>
  <tr>
    <th>m1</th>
    <th><tt>CATAGNCACGTGATCCTAGNCATGTTACCTGT</tt></th>
    <th></th>
  </tr>
  <tr>
    <th><i>A</i></th>
    <th><tt>01010001000010000100010000100000</tt></th>
    <th><tt>0x0422108a</tt></th>
  </tr>
  <tr>
    <th><i>C</i></th>
    <th><tt>10000010100000110000100000011000</tt></th>
    <th><tt>0x1810c141</tt></th>
  </tr>
  <tr>
    <th><i>G</i></th>
    <th><tt>00001000010100000010000100000010</tt></th>
    <th><tt>0x40840a10</tt></th>
  </tr>
  <tr>
    <th><i>T</i></th>
    <th><tt>00100000001001001000001011000101</tt></th>
    <th><tt>0xa3412404</tt></th>
    </tr>
  <tr>
    <th><i>A|C|G|T</i></th>
    <th><tt>11111011111111111110111111111111</tt></th>
    <th><tt>0xfff7ffdf</tt></th>
  </tr>
  </table>

We may set the above letter using <a href="https://software.intel.com/sites/landingpage/IntrinsicsGuide/">Intel Intrinsics</a>

<p>
  <tt>__m128i m1 = _mm_set_epi32(0xa3412404, 0x40840a10, 0x1810c141, 0x0422108a);</tt>
  </p>

<hr>

<b>Example 2</b>

Suppose we have two 32-symbol sequences (<tt>m1</tt> is shown above and <tt>m2</tt> is below).

<table>
  <tr>
    <th>m2</th>
    <th><tt>GCCTCAGTTTTCACTCTATCAATATGTAATAA</tt></th>
    <th></th>
  </tr>
  <tr>
    <th><i>A</i></th>
    <th><tt>00000100000010000100110100011011</tt></th>
    <th><tt>0xd8b21020</tt></th>
  </tr>
  <tr>
    <th><i>C</i></th>
    <th><tt>01101000000101010001000000000000</tt></th>
    <th><tt>0x0008a816</tt></th>
  </tr>
  <tr>
    <th><i>G</i></th>
    <th><tt>10000010000000000000000001000000</tt></th>
    <th><tt>0x02000041</tt></th>
  </tr>
  <tr>
    <th><i>T</i></th>
    <th><tt>00010001111000101010001010100100</tt></th>
    <th><tt>0x25454788</tt></th>
    </tr>
  <tr>
    <th><i>A|C|G|T</i></th>
    <th><tt>11111111111111111111111111111111</tt></th>
    <th><tt>0xffffffff</tt></th>
  </tr>
  </table>

We want to count the total number of symbols <tt>A</tt>, <tt>C</tt>, <tt>G</tt>, <tt>T</tt> that are at same positions for the both sequences. For this purcpose we perform a bitwise AND operation using <tt>_mm_and_si128</tt>, perform bitwise OR operation for <tt>A</tt>, <tt>C</tt>, <tt>G</tt>, <tt>T</tt> components and count the number of ones in the resulting 32-bit number using <tt>_mm_popcnt_u32</tt>.

<table>
  <tr>
    <th></th>
    <th></th>
    <th>A</th>
    <th>C</th>
    <th>G</th>
    <th>T</th>
  </tr>
  <tr>
    <th><tt>m1</tt></th>
    <th><tt>CATAGNCACGTGATCCTAGNCATGTTACCTGT</tt></th>
    <th><tt>0x0422108a</tt></th>
    <th><tt>0x1810c141</tt></th>
    <th><tt>0x40840a10</tt></th>
    <th><tt>0xa3412404</tt></th>
  </tr>
  <tr>
    <th><tt>m2</tt></th>
    <th><tt>GCCTCAGTTTTCACTCTATCAATATGTAATAA</tt></th>
    <th><tt>0xd8b21020</tt></th>
    <th><tt>0x0008a816</tt></th>
    <th><tt>0x02000041</tt></th>
    <th><tt>0x25454788</tt></th>
  </tr>
  <tr>
    <th><tt>m1 AND m2</tt></th>
    <th><tt>__________T_A__CTA___AT_T____T__</tt></th>
    <th><tt>0x00221000</tt></th>
    <th><tt>0x00008000</tt></th>
    <th><tt>0x00000000</tt></th>
    <th><tt>0x21410400</tt></th>
  </tr>
  <tr>
    <th colspan="6"><tt>0x00221000 OR 0x00008000 OR 0x00000000 OR 0x21410400 = 0x21639400</tt></th>
  </tr>
  <tr>
    <th colspan="6"><tt>_mm_popcnt_u32(0x21639400) = 9</tt></th>
  </tr>
</table>

<hr>



A code requires Intel's performance primitives to be used for sorting. It can be downloaded from <a>https://www.intel.com/content/www/us/en/developer/tools/oneapi/toolkits.html</a>

Human reference genome can be dowloaded from <a>https://www.ncbi.nlm.nih.gov/assembly/GCA_009914755.4</a>

The first step is to convert reference (FNA) and reads (FASTQ) files into binary files we need.

<h3>fna2bin</h3>
Convert the original FNA file into a binary file with index file (to store position of chunks)

<tt>fna2bin.exe C:\Genome\T2T\GCF_009914755.1_T2T-CHM13v2.0_genomic.fna C:\Genome\T2T\T2T_</tt>

<h4>Parameters</h4>

<ol>
  <li>Input FNA file</li>
  <li>Output folder + prefix</li>
</ol>

For the above example folder <tt>C:\Genome\T2T</tt> should exist. Two files will be created in the folder: <tt>T2T_data.bin</tt> and <tt>T2T_info.bin</tt>

File <tt>info</tt> contains the number of symbols in each chunk. The number of chunks can be found by dividing the file size by 4. It is assumed that FNA chunks contain no other symbols except <tt>A</tt>, <tt>C</tt>, <tt>G</tt>, <tt>T</tt>. Each symbols is coded by by two bits (<tt>A = 0</tt>, <tt>C = 1</tt>, <tt>G = 2</tt>, <tt>T = 3</tt>). Each chunk is rounded up to the nearest 32 symbols. So a chunk of 700 symbols will be rounded to <tt>22 x 32 = 704</tt> symbols and saved as <tt>22 x 8 = 176</tt> bytes in the binary file.

<h3>ref2m128</h3>
Convert the binary reference file (output of <b>fna2bin</b>) into a binary file supporting <tt>__m128i</tt> numbers.

<tt>ref2m128.exe C:\Genome\T2T\T2T_</tt>

<h4>Parameters</h4>

<ol>
  <li>Input/Output folder + prefix</li>
</ol>

Files <tt>T2T_data.bin</tt> and <tt>T2T_info.bin</tt> are in the folder <tt>C:\Temp2\Genome\T2T</tt>. File <tt>T2T_m128.bin</tt> is also created in the same folder (its size is twice big than <tt>T2T_data.bin</tt>).

<h3>fastq2bin</h3>
Convert a FASTQ file a binary file.

<tt>fastq2bin.exe D:\data\ERR263486_1.fastq C:\out\readsData.bin 100 10000000</tt>

<h4>Parameters</h4>

<ol>
  <li>Input FASTQ file</li>
  <li>Output binary file</li>
  <li>Length of reads</li>
  <li>Maximum number of reads</li>
</ol>

Only reads of a given length will be processed. A read is ignored if it containes any symbol except <tt>A</tt>, <tt>C</tt>, <tt>G</tt>, <tt>T</tt> (case insensitive). If a read is processed, then both read and its counterpart (flipped and <tt>A/T</tt>, <tt>C/G</tt> swapped) are written to the output file. Each symbol is coded by two bits (<tt>A = 0</tt>, <tt>C = 1</tt>, <tt>G = 2</tt>, <tt>T = 3</tt>). The number of bits for each read is rounded up to a byte. So if a read has 43 symbols, then we get 11 bytes in the output file.


<h3>bin2m128</h3>
Convert a binary FASTQ file (output of <b>fastq2bin</b>) to binary file supporting <tt>__m128i</tt> numbers. The output files will be used to perform alignment of reads.

<tt>bin2m128.exe C:\out\readsData.bin C:\out\readsData_m128.bin</tt>

<h4>Parameters</h4>

<ol>
  <li>Input binary file for reads</li>
  <li>Output binary file</li>
</ol>

<hr>
Once the binary files have been created, then for each read we find a list of candidate alignment positions with respect to the reference sequence, count the maximum number of matches attained for positions within the list, get a list of such poistions and count the number of its elements. There will be a lot of temporarily files created.It is better to have a temp folder to store all such files. Ideally, this folder should be on a fast drive (like SSD or NVMe) to speed of processing. The total size of files may be around 100~GB. As not all computers have large RAM memory, there is a parameter to control the number of chunks to be used.

Suppose we have already created the following folders:
<ol>
	<li><tt>inputRef</tt> (where reference files are, we use with the prefix for file names like <tt>C:\Genome\T2T\T2T_</tt>).</li>
	<li><tt>inputRead</tt> (path to binary reference file, both standard and m128)</li>
	<li><tt>tempFolder</tt> (all temporarily files will be stored there)</li>
	<li><tt>outputFolder</tt></li>
	<li><tt>outputPrefix</tt></li>
</ol>

Then we may use the followng sets of commands (described below in details).
<ol>
	<li><tt>createList.exe inputRef tempFolder seed</tt></li>
	<li><tt>sortList.exe tempFolder seed</tt></li>
<li><tt>searchPositions.exe inputRef tempFolder inputRead.bin tempFolder outputFolder/outputPrefix_pos.bin seed 0 999999999 4</tt></li>
<li><tt>countMatch.exe inputRef outputFolder/outputPrefix_pos.bin inputRead_m128.bin outputFolder/outputPrefix_match.bin</tt></li>
<li><tt>printMatch.exe outputFolder/outputPrefix_match.bin outputFolder/outputPrefix_stat.txt</tt></li>
</ol>

For example,
<ol>
<li>createList.exe D:\Data\Ref\T2T_ D:\Data\temp 111111111100011101100100100111010011100010100101000010100110000101111000000011</li>
<li>sortList.exe D:\Data\temp 111111111100011101100100100111010011100010100101000010100110000101111000000011</li>
<li>searchPositions.exe D:\Data\Ref\T2T_ D:\Data\temp D:\Data\Reads\err263486.bin D:\Data\temp D:\Data\out/S1_pos.bin 111111111100011101100100100111010011100010100101000010100110000101111000000011 0 99999 16</li>
<li>countMatch.exe D:\Data\Ref\T2T_ D:\Data\out/S1_pos.bin D:\Data\Reads\err263486_m128.bin D:\Data\out/S1_match.bin</li>
	<li>printMatch.exe D:\Data\out/S1_match.bin D:\Data\out/S1_stat.txt</li>
</ol>

<hr>

<h3>createList</h3>
Create a hash table for a given seed.

<tt>createList.exe C:\Genome\T2T\T2T_ D:\Temp\seed 1111111111000111011001001001110100111</tt>

<h4>Parameters</h4>

<ol>
  <li>Folder (+ prefix) for binary reference file</li>
  <li>Output folder</li>
  <li>Seed</li>
</ol>

A set of binary files is created (at most 65536 files) where parts of the hash table are stored. All files have <b>_orig.bin</b> suffix.



<h3>sortList</h3>
The hash table created by <b>createList</b> is sorted by "values". The sorting algorithm use Intel's Integrated Performance Primitives, so Intel software toolkit should be used to compile the code.

<tt>sortList.exe D:\Temp\seed 1111111111000111011001001001110100111</tt>

<h4>Parameters</h4>

<ol>
  <li>Input/output folder</li>
  <li>Seed</li>
</ol>

Hash pairs in the files are sorted by "values". For each original file with <b>_orig.bin</b> suffix a new file is created with <b>_sort.bin</b> suffix.

<h3>searchPositions</h3>
The hash table created by <b>createList</b> is sorted by "values". The sorting algorithm use Intel's Integrated Performance Primitives, so Intel software toolkit should be used to compile the code.

<tt>searchPositions.exe D:\Data\Ref\T2T_ D:\Temp\seed D:\Data\Reads\readData.bin D:\Temp\temp D:\out\readData_pos.bin 111011010100110101100100100101010110001 0 9999999 4</tt>

<h4>Parameters</h4>

<ol>
<li>Reference folder + prefix</li>
<li>Input folder (hash tables)</li>
<li>Input read binary file</li>
<li>Temp folder</li>
<li>Output file</li>
<li>Seed</li>
<li>First row to process</li>
<li>Last row to process</li>
<li>Number of chunks</li>
</ol>



<h3>countMatch</h3>
For a list of candidate positions we find the maximum number of matches as well as a list of positions where the maximum is attained.

<tt>countMatch.exe C:\Temp2\Genome\T2T\T2T_ D:\Gen2022\out\err263486_pos.bin D:\Gen2022\data\err263486_m128.bin D:\Gen2022\out\err263486_match.bin</tt>

<h4>Parameters</h4>

<ol>
<li>Reference folder + prefix</li>
<li>Input binary file (positions)</li>
<li>Input read file (m128)</li>
<li>Output binary file (matches)</li>
</ol>



