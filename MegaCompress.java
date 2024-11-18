// (setq fill-column 120)
//(setq fci-rule-column 120)
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Stack;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.FileReader;
import java.io.IOException;
import java.util.Collections;


public class MegaCompress
{
    public class hashTable
    {
	public class hashVectorElement
	{
	    public class preImageData
	    {
		/* If there are multiple strings that cause collisions we have one of these objects per string. The
		 * pre-image that the hash was computed from can be compared against str. */
		String str;
		/* Stores the start location of str in the string that str (the pre image) came from. Using this after
		 * the whole table is constructed we can replace the instances of str in the source string with an
		 * index. We can also use the size of this to see if we we gain anything from doing compression. */
		ArrayList<Integer> preImageOriginIndexes = new ArrayList<Integer>();

		preImageData(final String str, final int preImageOriginIndex)
		{
		    this.str = str;
		    preImageOriginIndexes.add(preImageOriginIndex);
		}

		/* For use when constructing a new table. Fills preImageOriginIndexes with multiple objects,
		 * representing multiple source locations. */
		preImageData(final String str, final ArrayList<Integer> preImageOriginIndexes)
		{
		    this.str = str;
		    this.preImageOriginIndexes = preImageOriginIndexes;
		}

		int getPreImageOriginIndexesSize()
		{
		    return preImageOriginIndexes.size();
		}
	    }
	
	    /* Stores any strings that produce hash. */
	    private ArrayList<preImageData> preImages = new ArrayList<preImageData>();

	    // ==========================================:( Member functions ):=========================================
	    public hashVectorElement(final String preImage, final int preImageOriginIndex)
	    {
		this.preImages.add(new preImageData(preImage, preImageOriginIndex));
	    }

	    /* For use when constructing a new table. Fills preImageOriginIndexes with multiple objects, representing
	     * multiple source locations. */
	    public hashVectorElement(final String preImage, final ArrayList<Integer> preImageOriginIndexes)
	    {
		this.preImages.add(new preImageData(preImage, preImageOriginIndexes));
	    }

	    /* This function adds a new element to preImages if there is a collision. Otherwise adds the index of the
	     * preImage in the original string to preImageOriginIndexes (this is assuming that preImageOriginIndex
	     * contains the right value.) */
	    public void potentialAddString(final String preImage, final int preImageOriginIndex)
	    {
		boolean found = false;
		for(int iter = 0; iter < preImages.size(); iter++)
		    {
			if(preImages.get(iter).str.equals(preImage))
			    {
				/* System.out.println("2nd copy of string added, with pre image origin of " +
				   preImageOriginIndex); */
				/* We've found a matching preImage, so all we need to do is update preImageOriginIndexes
				 * with the new pre image origin index. */
				preImages.get(iter).preImageOriginIndexes.add(preImageOriginIndex);
				found = true;
				break;
			    }
		    }
		if(!found)
		    {
			// System.out.println("Adding new pre image origin " + preImageOriginIndex);
			/* We didn't find the preImage in preImages. So we have a collision. We will update
			 * preImages. */
			preImages.add(new preImageData(preImage, preImageOriginIndex));
		    }
	    }

	    public void potentialAddString(final String preImage, final ArrayList<Integer> preImageOriginIndexes)
	    {
		boolean found = false;
		for(int iter = 0; iter < preImages.size(); iter++)
		    {
			if(preImages.get(iter).str.equals(preImage))
			    {
				System.out.println
				    ("Error: we shouldn't reach this point, because when we construct a new table " +
				     "we should only try to add the exact same string once, since preImages should " +
				     "only contain one copy of each string and we are copying from preImages.");
				System.exit(-1);
				// /* We've found a matching preImage, so all we
				//  * need to do is update preImageOriginIndexes
				//  * with the new pre image origin index. */
				// preImages.get(iter).preImageOriginIndexes
				//     = preImageOriginIndexes;
				// found = true;
				break;
			    }
		    }
		if(!found)
		    {
			// System.out.println("Adding pre image and pre image origins for new table pre image location");
			/* We didn't find the preImage in preImages. So we have a collision. We will update
			 * preImages. */
			preImages.add(new preImageData(preImage, preImageOriginIndexes));
		    }
	    }
	}
	
	/* Used to choose table size when expanding the table (the table should be a prime number in terms of it's size
	   because if the hash wraps around it will reduce clustering in terms of hits.) Primes contains all the primes
	   up to 1,000,000,000 with any primes that are more than 95% similar in magnitude removed (not including the
	   last prime.) The desired prime should be found by doubling the current table size and then using a binary
	   search to find the closest size that is larger than that. */
	private int [] primes =
	{2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 47, 53, 59, 61, 67, 71, 79, 83, 89, 97, 103, 107, 113, 127,
	 137, 139, 149, 151, 163, 167, 179, 181, 191, 193, 211, 223, 239, 241, 257, 263, 277, 281, 307, 311, 331, 337,
	 359, 367, 389, 397, 419, 421, 449, 457, 487, 491, 521, 523, 557, 563, 593, 599, 631, 641, 677, 683, 719, 727,
	 769, 773, 821, 823, 877, 881, 929, 937, 991, 997, 1051, 1061, 1117, 1123, 1187, 1193, 1259, 1277, 1361, 1367,
	 1439, 1447, 1531, 1543, 1627, 1637, 1733, 1741, 1847, 1861, 1973, 1979, 2087, 2089, 2203, 2207, 2333, 2339,
	 2467, 2473, 2609, 2617, 2767, 2777, 2927, 2939, 3109, 3119, 3299, 3301, 3491, 3499, 3691, 3697, 3907, 3911,
	 4127, 4129, 4349, 4357, 4591, 4597, 4861, 4871, 5147, 5153, 5431, 5437, 5737, 5741, 6047, 6053, 6373, 6379,
	 6719, 6733, 7103, 7109, 7487, 7489, 7901, 7907, 8329, 8353, 8803, 8807, 9277, 9281, 9781, 9787, 10303, 10313,
	 10859, 10861, 11437, 11443, 12049, 12071, 12713, 12721, 13397, 13399, 14107, 14143, 14891, 14897, 15683, 15727,
	 16561, 16567, 17443, 17449, 18371, 18379, 19373, 19379, 20399, 20407, 21487, 21491, 22637, 22639, 23831, 23833,
	 25097, 25111, 26437, 26449, 27847, 27851, 29327, 29333, 30881, 30893, 32531, 32533, 34253, 34259, 36067, 36073,
	 37987, 37991, 40009, 40013, 42131, 42139, 44357, 44371, 46723, 46727, 49193, 49199, 51797, 51803, 54539, 54541,
	 57413, 57427, 60457, 60493, 63689, 63691, 67049, 67057, 70589, 70607, 74353, 74357, 78277, 78283, 82421, 82457,
	 86813, 86837, 91411, 91423, 96259, 96263, 101333, 101341, 106681, 106693, 112327, 112331, 118247, 118249,
	 124477, 124489, 131059, 131063, 137983, 137993, 145259, 145267, 152939, 152941, 160997, 161009, 169489, 169493,
	 178417, 178439, 187843, 187861, 197753, 197759, 208189, 208207, 219169, 219187, 230729, 230743, 242911, 242923,
	 255709, 255713, 269177, 269179, 283369, 283397, 298327, 298339, 314059, 314063, 330607, 330611, 348017, 348031,
	 366383, 366397, 385709, 385739, 406067, 406073, 427447, 427451, 449951, 449959, 473647, 473659, 498599, 498611,
	 524857, 524863, 552491, 552493, 581573, 581597, 612217, 612223, 644447, 644489, 678413, 678421, 714139, 714143,
	 751739, 751747, 791317, 791321, 832973, 832987, 876833, 876851, 923017, 923023, 971639, 971651, 1022797,
	 1022821, 1076657, 1076671, 1133357, 1133359, 1193011, 1193021, 1255829, 1255831, 1321939, 1321951, 1391549,
	 1391557, 1464809, 1464811, 1541921, 1541923, 1623077, 1623091, 1708523, 1708529, 1798457, 1798469, 1893131,
	 1893163, 1992817, 1992839, 2097727, 2097743, 2208203, 2208229, 2324453, 2324459, 2446811, 2446813, 2575627,
	 2575633, 2711197, 2711213, 2853911, 2853913, 3004123, 3004151, 3162277, 3162283, 3328723, 3328729, 3503933,
	 3503939, 3688357, 3688361, 3882523, 3882539, 4086889, 4086911, 4302013, 4302017, 4528487, 4528507, 4766863,
	 4766887, 5017811, 5017819, 5281919, 5281921, 5559923, 5559929, 5852569, 5852591, 6160633, 6160639, 6484889,
	 6484897, 6826271, 6826289, 7185571, 7185623, 7563821, 7563833, 7961939, 7961951, 8381003, 8381047, 8822159,
	 8822161, 9286499, 9286507, 9775307, 9775313, 10289819, 10289821, 10831391, 10831423, 11401549, 11401561,
	 12001667, 12001697, 12633367, 12633373, 13298323, 13298333, 13998251, 13998301, 14735081, 14735093, 15510641,
	 15510643, 16326997, 16327007, 17186341, 17186347, 18090899, 18090901, 19043081, 19043159, 20045447, 20045449,
	 21100489, 21100517, 22211099, 22211131, 23380139, 23380163, 24610699, 24610711, 25906057, 25906061, 27269563,
	 27269573, 28704883, 28704889, 30215677, 30215699, 31805999, 31806001, 33480023, 33480043, 35242183, 35242213,
	 37097083, 37097129, 39049613, 39049649, 41104897, 41104907, 43268353, 43268369, 45545653, 45545657, 47942801,
	 47942803, 50466167, 50466193, 53122327, 53122339, 55918259, 55918277, 58861349, 58861379, 61959367, 61959379,
	 65220401, 65220437, 68653097, 68653103, 72266431, 72266449, 76070021, 76070023, 80073731, 80073737, 84288173,
	 84288229, 88724453, 88724459, 93394183, 93394211, 98309707, 98309719, 103483943, 103483951, 108930517,
	 108930533, 114663719, 114663737, 120698671, 120698681, 127051273, 127051277, 133738211, 133738219, 140777083,
	 140777093, 148186421, 148186429, 155985733, 155985749, 164195621, 164195641, 172837519, 172837541, 181934273,
	 181934287, 191509793, 191509853, 201589327, 201589331, 212199301, 212199307, 223367747, 223367761, 235123961,
	 235123969, 247498949, 247498969, 260525231, 260525249, 274237121, 274237123, 288670673, 288670687, 303863909,
	 303863927, 319856827, 319856839, 336691427, 336691451, 354412061, 354412081, 373065379, 373065383, 392700433,
	 392700439, 413368897, 413368913, 435125189, 435125221, 458026549, 458026573, 482133257, 482133259, 507508699,
	 507508709, 534219733, 534219737, 562336571, 562336603, 591933269, 591933283, 623087669, 623087681, 655881773,
	 655881781, 690401897, 690401923, 726738877, 726738893, 764988319, 764988353, 805250903, 805250923, 847632571,
	 847632581, 892244833, 892244839, 939205123, 939205147, 988637033, 988637053, 999999937};
	private int initialTableSize = 113;
	private int currentTableLoadingFactor = 0;
	/* We should expand the table if currentTableLoadingFactor exceeds this value. */
	private double expandTableThreshold = 0.7;
	private hashVectorElement [] table = new hashVectorElement [initialTableSize];

	private int hash(final String str, final int tblSize)
	{
	    final int prime = 37;
	    int hash = 0;

	    for(char c: str.toCharArray())
		{
		    hash = hash * prime + c;
		}
	    // Unsigned integer don't exist in Java >:(.
	    if(hash < 0)
		{
		    hash = -hash;
		}
	    // System.out.println("hash = " + (int)(hash % tblSize) + ", for str \"" + str + "\"");
	    return (int)(hash % tblSize);
	}

	public void add(final String preImage, final int preImageOriginIndex)
	{
	    final int preImageHash = hash(preImage, table.length);
	    
	    if(table[preImageHash] == null)
		{
		    // System.out.println("No collision.");
		    /* We're adding a new item to table, we must inc our load factor variable. */
		    currentTableLoadingFactor++;
		    table[preImageHash] = new hashVectorElement
			(preImage, preImageOriginIndex);
		}
	    else
		{
		    // System.out.println("PreImageHash collision.");
		    table[preImageHash].potentialAddString(preImage, preImageOriginIndex);
		}

	    if((double)currentTableLoadingFactor / table.length > expandTableThreshold)
		{
		    /* Our table has read a load factor of expandTableThreshold and so we will expand it. */
		    final int newDesiredSizeApproximation = table.length * 2;
		    
		    if(newDesiredSizeApproximation > primes[primes.length -1])
			{
			    System.out.println
				("Error: table in object of type hashTable has exceeded it's maximum loading factor " +
				 "and needs resizing, however when trying to resize it has been found that the " +
				 "estimation for the new size is greater than the value of the largest entry (" +
				 primes[primes.length -1] + ") in the pre-computed primes table.");
			    System.exit(-1);
			}
		    else
			{
			    int lowPos = 0;
			    int midPos = 0;
			    int highPos = primes.length -1;
			    
			    while(lowPos <= highPos)
				{
				    midPos = lowPos + (highPos - lowPos) / 2;

				    if(primes[midPos] == newDesiredSizeApproximation)
					{
					    // Found;
					    break;
					}
				    else if(primes[midPos] < newDesiredSizeApproximation)
					{
					    // Search to the right.
					    lowPos = midPos + 1;
					}
				    else
					{
					    // Search to the left.
					    highPos = midPos -1;
					}
				}

			    int newTableSize = 0;
			    
			    if(primes[midPos] != newDesiredSizeApproximation)
				{
				    /* We already checked above if newDesiredSizeApproximation was larger than the
				       maximum element in primes. So if it's not equal there should be at least one more
				       element after midPos. */
				    newTableSize = primes[midPos + 1];
				}
			    else
				{
				    // An exact match! すごいですね!
				    newTableSize = newDesiredSizeApproximation;
				}

			    currentTableLoadingFactor = 0;
			    hashVectorElement [] newTable = new hashVectorElement [newTableSize];

			    for(int iter = 0; iter < table.length; iter++)
				{
				    if(table[iter] != null)
					{
					    /* This location has something in it, so we will get all of the preImages
					     * and use them to calculate new hashes to construct the new table. We also
					     * need to use the preImages to find the indexes into stringVector for the
					     * preImages. */
					    ArrayList<hashVectorElement.preImageData> preImages = table[iter].preImages;

					    for(int preImageIter = 0; preImageIter < preImages.size(); preImageIter++)
						{
						    final String preImageLocal = preImages.get(preImageIter).str;
						    final ArrayList<Integer> preImageOriginIndexes =
							preImages.get(preImageIter).preImageOriginIndexes;
						    // Compute new hash.
						    final int preImageHashLocal = hash(preImageLocal, newTable.length);
						    if(newTable[preImageHashLocal] == null)
							{
							    currentTableLoadingFactor++;
							    newTable[preImageHashLocal] = new hashVectorElement
								(preImageLocal, preImageOriginIndexes);
							}
						    else
							{
							    newTable[preImageHashLocal].potentialAddString
								(preImageLocal, preImageOriginIndexes);
							}
						}
					}
				}

			    table = newTable;
			}
		}
	}
	int size()
	{
	    return table.length;
	}

	// Can return a null object.
	hashVectorElement get(final int index)
	{
	    if(index < 0 || index > table.length -1)
		{
		    System.out.println
			("Error invalid index (" + index +
			 ") passed to hashTable.get(), should be in the " +
			 "range [0, " + table.length + ")." );
		    System.exit(-1);
		}

	    return table[index];
	}


	/* Returns a sorted (and expanded) vector of any elements in the hash table that are non-null. The vector is
	 * expanded in the sense that there is an entry for every pre-image. That is for every non-null element in table
	 * there will be an entry for every pre-image entry in preImages that it has. The returned vector will be sorted
	 * (in descending order) by the preImageOriginIndexes.size() value for each entry in preImages, in descending
	 * order. */
	ArrayList<hashVectorElement.preImageData> getElementsSortedByPreImageOccurrenceCount()
	{
	    ArrayList<hashVectorElement.preImageData> preImagesRet =  new ArrayList<hashVectorElement.preImageData>();
	    for(hashVectorElement hVE: table)
		{
		    if(hVE != null)
			{
			    for(hashVectorElement.preImageData pIs: hVE.preImages)
				{
				    preImagesRet.add(pIs);
				}
			}
		}
	    preImagesRet.sort(Comparator.comparing
			      (hashVectorElement.preImageData::getPreImageOriginIndexesSize).reversed());
	    return preImagesRet;
	}
    }

    public static String minEncodedSizeWithRLECompression(String sequence, int maxUnitSize)
    {
	class dynamicProgrammingInfo
	{
	    /* Compression size of stream up to the current point. This is basically used to calculate if we get any
	     * advantage from the current compression sequence that has been calculated. */
	    int sizeToPoint;
	    /* Index of point that leads to this point in the compression stream. */
	    int lastPoint = 0;
	    // Size of the units that are repeated.
	    int runLengthUnitSize = 0;;
	    /* Number of units that are to be repeated. If the repeatCount is not > 1, then there is no compression for
	     * this point. */
	    int repeatCount = 0;

	    public String toString()
	    {
		return sizeToPoint + ", " + lastPoint +
		    ", " + runLengthUnitSize + ", " + repeatCount;
	    }
	}
	
        int seqLen = sequence.length();
        dynamicProgrammingInfo [] dp = new dynamicProgrammingInfo [seqLen];

        for(int i = 0; i < seqLen; i++)
	    {
		dp[i] = new dynamicProgrammingInfo();
		/* Initialize dp array positions to their index (as the dp array stores the current best compression
		 * length found up to that point in the character stream to be compressed and so the current index is
		 * the worse case (no compression.)) */
		dp[i].sizeToPoint = i + 1;
		/* Make sure that we are pointing to the last point in the stream. If there is no compression for this
		   point we still want to be able to get back to the point before it when we work backwards through the
		   calculated compression sequence to the start. Note that we want dp[0].lastPoint to be 0 as a special
		   case. */
		dp[i].lastPoint = i;
	    }
        
        // Dynamic Programming to calculate minimum encoding size
        for(int i = 0; i < seqLen; i++)
	    {
		// Try different run lengths from 1 to maxUnitSize
		for(int runLength = 1; runLength <= maxUnitSize; runLength++)
		    {
			if(i + runLength < seqLen)
			    {
				String unit = sequence.substring(i, i + runLength);
				int repeatCount = countRepeats(sequence, unit, i);
				/* We don't need ...repeatCount).length() as seen above because we only need one thing
				 * for the start of the run length and another for the number (since we are storing in
				 * binary.) That's 2 + the length of the unit that will be repeated. */
				int encodingSize = 2 + unit.length();
				int uncompressedIndex = i + runLength * repeatCount -1;

				/* Update dp array and track decision. In this if statment we first check if we will do
				 * an out of bounds index and if not we check the best compression size from the current
				 * position dp[i] plus the newly calculated encoding size against the current best
				 * compression size for the position at the end of the sequence we just calculated the
				 * compression for dp[uncompressedIndex] If the compression encoding is smaller we
				 * store that this is the best compression we have found upto this point and also where
				 * the compression starts (we will probably need to store the size of the unit that we
				 * are repeating too (runLength). */
				if(repeatCount > 1 &&
				   uncompressedIndex <= seqLen &&
				   dp[i].sizeToPoint -1 + encodingSize <
				   dp[uncompressedIndex].sizeToPoint)
				    {
					dp[uncompressedIndex].sizeToPoint = dp[i].sizeToPoint -1 + encodingSize;
					/* Make sure we know what point was before this one in the sequence so we can
					   reach the start of the sequence later. */
					dp[uncompressedIndex].lastPoint = i;
					dp[uncompressedIndex].runLengthUnitSize = runLength;
					dp[uncompressedIndex].repeatCount = repeatCount;
				    }
			    }
		    }
	    }

	Stack<dynamicProgrammingInfo> compressionSequence = new Stack<>();
	/* Back track through dp array following the compression sequence that was calculated (while saving the sequence
	 * in the correct order.) Note that dp[0].lastPoint should be 0. */
	for(int dpIndex = seqLen -1; dpIndex > -1; )
	    {
		compressionSequence.push(dp[dpIndex]);
		dpIndex = dp[dpIndex].lastPoint -1;
	    }

	// Now, reconstruct the compressed sequence using the decision array
        StringBuilder compressedSequence = new StringBuilder();
	int atPointInSequence = 0;

        while(!compressionSequence.isEmpty())
	    {
		/* If the repeat count is less than 2 there is no compression at this point. */
		if(compressionSequence.peek().repeatCount > 1)
		    {
			int atPointInSequenceLocal = atPointInSequence;
			// Get unit that will be repeated.
			StringBuilder rLUnit = new StringBuilder();
			for(int iter = 0; iter < compressionSequence.peek().runLengthUnitSize; iter++)
			    {
				rLUnit.append(sequence.charAt(atPointInSequenceLocal));
				atPointInSequenceLocal++;
			    }
			/* Add RL encoded sequence to compressedSequence. */
			compressedSequence.append("ラ" + compressionSequence.peek().runLengthUnitSize);
			compressedSequence.append(rLUnit);
			compressedSequence.append(Integer.toString(compressionSequence.peek().repeatCount));
			/* Update atPointInSequence (which is needed when adding chars without compression (see else
			 * below.) */
			atPointInSequence += compressionSequence.peek().runLengthUnitSize *
			    compressionSequence.peek().repeatCount;
		    }
		else
		    {
			/* Push back the current char and advance one point in the sequence. */
			compressedSequence.append(sequence.charAt(atPointInSequence));
			atPointInSequence++;
		    }
		compressionSequence.pop();
	    }

        return compressedSequence.toString();
    }

    // Helper method to count consecutive repeats of a substring from a given start
    private static int countRepeats(final String sequence, final String unit, final int start)
    {
        int count = 0;
        int i = start;
        while (i < sequence.length() && sequence.startsWith(unit, i))
	    {
		count++;
		i += unit.length();
	    }
        return count;
    }


    public static String uncompressRLE(final String rLECompressedInput)
    {
	StringBuilder uncompressedRet = new StringBuilder();
	for(int iter = 0; iter < rLECompressedInput.length(); iter++)
	    {
		if(rLECompressedInput.charAt(iter) == 'ラ')
		    {
			StringBuilder unit = new StringBuilder();
			iter++;

			// Get unit size.
			int unitSize = 0;
			while(iter < rLECompressedInput.length() && Character.isDigit(rLECompressedInput.charAt(iter)))
			    {
				unitSize = unitSize * 10 + (rLECompressedInput.charAt(iter) - '0');
				iter++;
			    }
			
			for(int uIter = 0; uIter < unitSize; uIter++)
			    {
				unit.append(rLECompressedInput.charAt(iter + uIter));
			    }
			iter += unitSize;

			// Get run length size.
			int rL = 0;
			while(iter < rLECompressedInput.length() && Character.isDigit(rLECompressedInput.charAt(iter)))
			    {
				rL = rL * 10 + (rLECompressedInput.charAt(iter) - '0');
				iter++;
			    }
			iter--;

			for(int rLIter = 0; rLIter < rL; rLIter++)
			    {
				uncompressedRet.append(unit);
			    }
		    }
		else
		    {
			uncompressedRet.append(rLECompressedInput.charAt(iter));
		    }
	    }

	return uncompressedRet.toString();
    }

    /* Minimum size that a reference takes up ("JNN", where N is is '0' to '9' and J is a Japanese character.).
       Since maxRefIndexNumber is 2x maxRefIndexNumberComponentPerChar, it should fit in two characters. */
    static final int referenceSize = 3;
    static final int maxRefIndexNumberComponentPerChar = 128;
    static final int maxRefIndexNumber = maxRefIndexNumberComponentPerChar * maxRefIndexNumberComponentPerChar;
    
    public String doDictionaryCompressionOnRLERepresentation
	(String rLECompressedInput, final char compressionSignifier, final char compressionDictEndSignifier,
	 final int subStrLen)
    {
	String dictTable = "";
	String currentSubStr = "";
	hashTable localHashTable = new hashTable();

	if(subStrLen < referenceSize + 1)
	    {
		System.out.println
		    ("Error: value of subStrLen (" + subStrLen + ") passed to " +
		     "doDictionaryCompressionOnRLERepresentation() is less than the minimum maximum reference size (" +
		     referenceSize + ").");
		System.exit(-1);
	    }

	// ======================================================== Build hash table! ==================================
	for(int iter = 0; iter < rLECompressedInput.length() - subStrLen; iter++)
	    {
		final String preImage = rLECompressedInput.substring(iter, iter + subStrLen);
		localHashTable.add(preImage, iter);
	    }

	// Should be rLECompressedInput.length() in size.
        char [] dictCompressedInput		= rLECompressedInput.toCharArray();
	/* Map of which locations we've altered so far (we need this to check for possible overlapping repeated strings,
	 * as we can't write over locations that have already been altered (in the general case.)) */
	final char unAlteredInputMapChar	= ' ';
	final char alteredInputMapChar	 	= 'A';
	char [] dictCompressedInputAlteredMap 	= new char [rLECompressedInput.length()];
        Arrays.fill(dictCompressedInputAlteredMap, unAlteredInputMapChar); // Fill the array with ' '.

	ArrayList<String> stringVector = new ArrayList<String>();

	int tmpSafePreImageOriginIndexesSize = 0;

	// ======================================================== Iterate through all valid entries in the hash table.
	for(hashTable.hashVectorElement.preImageData pID: localHashTable.getElementsSortedByPreImageOccurrenceCount())
	    {		
		if(pID.getPreImageOriginIndexesSize() > 1)
		    {
			/* This variable should be updated when (and only when) we first add the current
			 * pre-image to the stringVector. */
			int preImageIndexInStringVector = 0;
			
			final int preImgLen = pID.str.length();


			/* Contains indexes for the current pre-image that don't conflict with what has already been
			 * written to dictCompressedInputAlteredMap. */
			ArrayList<Integer> safePreImageOriginIndexes = new ArrayList<Integer>();

			/* ======= Iterate through all origins (in the source string) of this pre-image and save indexes
			 * that don't point to locations that have already been written to. If there are more than 2,
			 * then we may benefit from replacing the pre-images with the indexes (assuming that they don't
			 * overlap with each other.)  */
			for(int preImageOriginLocationsIter = 0;
			    preImageOriginLocationsIter < pID.preImageOriginIndexes.size();
			    preImageOriginLocationsIter++)
			    {
				final int preImageSourceIndex =  pID.preImageOriginIndexes.get
				    (preImageOriginLocationsIter);
						
				/* Check if we've already updated this location. If so then we will skip updating it, so
				 * as not to mess up what was already done. */
				boolean notAlreadyDiddled = true;
				for(int preImgStrIter = 0; preImgStrIter < preImgLen; preImgStrIter++)
				    {
					if(dictCompressedInputAlteredMap[preImageSourceIndex + preImgStrIter] !=
					   unAlteredInputMapChar)
					    {
						/* Don't add a reference because it would "mask" a previous
						 * reference. */
						notAlreadyDiddled = false;
						break;
					    }
				    }
				if(notAlreadyDiddled)
				    {
					safePreImageOriginIndexes.add(preImageSourceIndex);
				    }
			    }


			/* ======= Check for overlapping pre-image origins for the current pre-image origin indexes and
			 * length and remove overlapping pre-images. This won't always find an optimal arrangement,
			 * however it will remove potential errors. */
			/* NOTE: SINCE EACH PRE-IMAGE WAS CHECKED FOR IN ORDER safePreImageOriginIndexes SHOULD BE IN
			 * ORDER. IF IT IS NOT THE THE FOLLOWING CODE WILL NOT BE CORRECT! SEE
			 * localHashTable.add(preImage, iter); ELSEWHERE (ASSUMING IT'S NOT BEEN CHANGED.) */
			if(safePreImageOriginIndexes.size() > 2)
			    {
				ArrayList<Integer> fullySafePreImageOriginIndexes = new ArrayList<Integer>();
				/* We always add the first pre-image as we consider the conflicting image to be the next
				 * one (if any.) */
				fullySafePreImageOriginIndexes.add(safePreImageOriginIndexes.get(0));
				    
				for(int preImageOriginLocationsIter = 0;
				    preImageOriginLocationsIter < safePreImageOriginIndexes.size() -1;
				    preImageOriginLocationsIter++)
				    {
					/* Since (as stated above) we have assumed that safePreImageOriginIndexes is in
					 * order we only need to check if the end of the pre-image instance starting at
					 * preImageOriginLocationsIter is equal to or more than the start of the
					 * pre-image instance starting at preImageOriginLocationsIter + 1. */
					if(safePreImageOriginIndexes.get(preImageOriginLocationsIter) + preImgLen <
					   safePreImageOriginIndexes.get(preImageOriginLocationsIter + 1))
					    {
						fullySafePreImageOriginIndexes.add
						    (safePreImageOriginIndexes.get(preImageOriginLocationsIter + 1));
					    }
				    }
				    
				safePreImageOriginIndexes = fullySafePreImageOriginIndexes;
			    }

			/* ======= Replace pre-image instances pointed to by safePreImageOriginIndexes with references
			 * to the pre-image. By the following if statement we know if we will benefit from replacing the
			 * pre-image instances with references. We need +1 below to account for compressionSignifier. */
			if((referenceSize * safePreImageOriginIndexes.size() + subStrLen + 1) <
			   (subStrLen * safePreImageOriginIndexes.size()))
			    {
				tmpSafePreImageOriginIndexesSize += safePreImageOriginIndexes.size();
				/* Add pre-image to stringVector! */
				preImageIndexInStringVector = stringVector.size();
				stringVector.add(pID.str);

				/* Format the index of the pre image (in the source string) as a
				   referenceSize-character -1 wide string with leading zeros. */
				if(preImageIndexInStringVector > maxRefIndexNumber)
				    {
					System.out.println
					    ("Error: preImageIndexInStringVector (" + preImageIndexInStringVector +
					     ") " + "found to be larger than maxRefIndexNumber (" + maxRefIndexNumber +
					     ").");
					System.exit(-1);
				    }
				
				ArrayList<Character> stringVectorIndexStr = new ArrayList<>();
				/* Convert the string vec index number to base-128 and then add 128.  */
				if(preImageIndexInStringVector > 0)
				    {
					while(preImageIndexInStringVector > 0)
					    {
						// Get remainder.
						int rem = preImageIndexInStringVector % maxRefIndexNumberComponentPerChar;
						// Map to range (127, 255].
						char c = (char)(rem + maxRefIndexNumberComponentPerChar);
						stringVectorIndexStr.add(c);
						// Get the quotient.
						preImageIndexInStringVector /= maxRefIndexNumberComponentPerChar;
					    }
				    }
				else
				    {
					for(int iter = 0; iter < referenceSize - 1; iter++)
					    {
						stringVectorIndexStr.add((char)maxRefIndexNumberComponentPerChar);
					    }
				    }
				if(stringVectorIndexStr.size() < referenceSize - 1)
				    {
					// Add 0 padding.
					for(int iter = 0; (referenceSize - 1) - stringVectorIndexStr.size() > 0; iter++)
					    {
						stringVectorIndexStr.add((char)maxRefIndexNumberComponentPerChar);
					    }
				    }
				/* StringVectorIndexStr is populated by the least significant digit first. So it
				   will have to be reversed. */
				Collections.reverse(stringVectorIndexStr);
				
				/* ======= Iterate through all origins safe (not already written over) locations (in the
				 * source string) for the current pre-image and insert references. */
				for(int preImgSrcIndex: safePreImageOriginIndexes)
				    {
					/* compressionSignifier will indicate that a index follows. */
					dictCompressedInput[preImgSrcIndex]			= compressionSignifier;
					dictCompressedInputAlteredMap[preImgSrcIndex] 		= stringVectorIndexStr.get(0);
					dictCompressedInput[preImgSrcIndex + 1]			= stringVectorIndexStr.get(0);
					dictCompressedInputAlteredMap[preImgSrcIndex + 1] 	= stringVectorIndexStr.get(1);
					dictCompressedInput[preImgSrcIndex + 2] 		= stringVectorIndexStr.get(1);							
					for(int preImgStrIter = referenceSize; preImgStrIter < preImgLen;
					    preImgStrIter++)
					    {
						dictCompressedInput[preImgSrcIndex + preImgStrIter]	= '`';
						dictCompressedInputAlteredMap
						    [preImgSrcIndex  + preImgStrIter]	= alteredInputMapChar;
					    }
				    }
			    }
		    }
		else
		    {
			/* Since localHashTable.getElementsSortedByPreImageOccurrenceCount() returns a sorted list, if
			   we encounter any pID where pID.getPreImageOriginIndexesSize() = 1, we know that any other
			   pIDs that remain to be checked will also have getPreImageOriginIndexesSize() return 1. */
			break;
		    }
	    }

	if(stringVector.size() > 0)
	    {
		StringBuilder dictCompressedInputProper = new StringBuilder();
		// Add string data that will be referenced.
		for(int iter = 0; iter < stringVector.size(); iter++)
		    {
			dictCompressedInputProper.append(stringVector.get(iter) + compressionSignifier);
		    }
		/* Add special end character to denote the start of the data and references. compressionDictEndSignifier
		   represents the start of the main data and the end of the list of strings. */
		dictCompressedInputProper.setCharAt(dictCompressedInputProper.length() -1, compressionDictEndSignifier);
		// Create compacted representation!
		for(int iter = 0; iter < dictCompressedInput.length; iter++)
		    {
			if(dictCompressedInput[iter] != '`')
			    {
				dictCompressedInputProper.append(dictCompressedInput[iter]);
			    }
		    }

		return dictCompressedInputProper.toString();
	    }
	else
	    {
		return rLECompressedInput;
	    }
    }


    public static String uncompressDictionaryCompression
	(String dictionaryCompressedStream, final char compressionSignifier, final char compressionDictEndSignifier)
    {
	ArrayList<String> stringVector = new ArrayList<String>();

	// Build stringVector...
	StringBuilder string = new StringBuilder();
	int iter = 0;
	boolean foundEndOfDict = false;

	/* Note here that we assume that the dictionary section contains at least one entry before
	 * compressionDictEndSignifier. A well formed input should have this property anyway. */
	for( ; iter < dictionaryCompressedStream.length() ; iter++)
	    {
		if(dictionaryCompressedStream.charAt(iter) == compressionDictEndSignifier)
		    {
			foundEndOfDict = true;
			break;
		    }
		
		if(dictionaryCompressedStream.charAt(iter) == compressionSignifier)
		    {
			stringVector.add(string.toString());
			string = new StringBuilder();
		    }
		else
		    {
			string.append(dictionaryCompressedStream.charAt(iter));
		    }
	    }
	// Add last string.
	stringVector.add(string.toString());

	if(foundEndOfDict)
	    {
		// Move past compressionDictEndSignifier.
		iter++;
		string = new StringBuilder();
		/* Add some error message here and exit if we've reached the end of the steam? */
		
		for( ; iter < dictionaryCompressedStream.length(); iter++)
		    {
			if(dictionaryCompressedStream.charAt(iter) == compressionSignifier)
			    {
				int strVecIndex = 0;
				// Move past compressionSignifier.
				iter++;
				for(int indexLen = 0; indexLen < referenceSize - 1; indexLen++, iter++)
				    {
					if(iter == dictionaryCompressedStream.length() -1)
					    {
						// Put a better error message here.
						System.out.println
						    ("Error: reached end of stream before expected.");
						System.exit(-1);
					    }
					
					// Get the character containing the number.
					char c = dictionaryCompressedStream.charAt(iter);
					/* Extract the base-maxRefIndexNumberComponentPerChar val
					   (subtract maxRefIndexNumberComponentPerChar to get the original remainder,
					   since it had that much added to it.) */
					int val = c - maxRefIndexNumberComponentPerChar;
					/* Shift the current number by maxRefIndexNumberComponentPerChar (move it one
					   "digit" to the left in base-maxRefIndexNumberComponentPerChar) */
					strVecIndex = strVecIndex * maxRefIndexNumberComponentPerChar + val;
				    }
				/* We've moved to the character after the end of the dictionary index number, however
				 * the outer for loop will also increment iter. So we must step back 1 here to account
				 * for that. */
				iter--;
				
				try 
				    {
					string.append(stringVector.get(strVecIndex));
				    }
				catch (Exception e) 
				    {
					System.out.println(dictionaryCompressedStream.substring(iter - 20, iter + 20));
					System.out.println(stringVector.size());
					System.out.println(strVecIndex);				
					// Log the exception, perform some cleanup, etc.
					throw e; // Re-throw the caught exception
				    }				
			    }
			else
			    {
				string.append(dictionaryCompressedStream.charAt(iter));
			    }
		    }

		return string.toString();
	    }
	else
	    {
		return dictionaryCompressedStream;
	    }
    }

    
    public static void main(String[] args)
    {
	MegaCompress MegaCompressInstance = new MegaCompress();

	// Read in sequencePath file into sequence.
	String sequence = "";
	{
	    String sequencePath = "preimageText.txt";
	    StringBuilder sBSequence = new StringBuilder();
	    try(BufferedReader br = new BufferedReader(new FileReader(sequencePath)))
		{
		    String line;
		    while((line = br.readLine()) != null)
			{
			    sBSequence.append(line).append(System.lineSeparator()); // Append line breaks
			}
		}
	    catch(IOException e)
		{
		    System.out.println("Error: failed to read file \"" + sequencePath + "\"");
		    e.printStackTrace();
		    System.exit(-1);
		}
	    sequence = sBSequence.toString();
	}
	
	
        int maxUnitSize = 20; // Maximum length of the unit that can be repeated
        String rLECompressed = minEncodedSizeWithRLECompression(sequence, maxUnitSize);

	    
	char [] dictCompressionSignifiers		=
	    {'あ', 'い', 'う', 'え', 'お', 'か', 'き', 'く', 'け', 'こ', 'さ', 'し', 'す', 'せ', 'そ', 'た', 'ち', 'つ',
	     'て', 'と', 'な', 'に', 'ぬ', 'ね', 'の', 'は', 'ひ', 'ふ', 'へ', 'ほ', 'ま', 'み', 'む', 'め', 'も', 'や',
	     'ゆ', 'よ', 'ら', 'り', 'る', 'れ', 'ろ', 'わ', 'ゐ', 'ゑ'};
	char [] dictCompressionDictEndSignifiers	=
	    {'ア', 'イ', 'ウ', 'エ', 'オ', 'カ', 'キ', 'ク', 'ケ', 'コ', 'サ', 'シ', 'ス', 'セ', 'ソ', 'タ', 'チ', 'ツ',
	     'テ', 'ト', 'ナ', 'ニ', 'ヌ', 'ネ', 'ノ', 'ハ', 'ヒ', 'フ', 'ヘ', 'ホ', 'マ', 'ミ', 'ム', 'メ', 'モ', 'ヤ',
	     'ユ', 'ヨ', 'リ', 'ル', 'レ', 'ロ', 'ワ', 'ヰ', 'ヱ', 'ヲ'};
	
	String dictAndRLECompressed = rLECompressed;
	for(int iter = dictCompressionSignifiers.length -1; iter >= 0; iter--)
	    {
		/* We must use MaxReferenceSize + 1, since the sub-string should be larger than the max reference size
		 * to get any gain.  */
		dictAndRLECompressed = MegaCompressInstance.doDictionaryCompressionOnRLERepresentation
		    (dictAndRLECompressed, dictCompressionSignifiers[iter],
		     dictCompressionDictEndSignifiers[iter], referenceSize + 1 + iter);
	    }
	
	String rLEUncompressed = uncompressRLE(rLECompressed);
	String dictAndRLEUnCompressed = dictAndRLECompressed;

	for(int iter = 0; iter < dictCompressionSignifiers.length; iter++)
	    {
		dictAndRLEUnCompressed = MegaCompressInstance.uncompressDictionaryCompression
		     (dictAndRLEUnCompressed, dictCompressionSignifiers[iter], dictCompressionDictEndSignifiers[iter]);
	    }

	dictAndRLEUnCompressed = uncompressRLE(dictAndRLEUnCompressed);

	// System.out.println("=========================RLECompressed = \n" + rLECompressed);
	// System.out.println("\n=========================RLEUnCompressed = \n" + rLEUncompressed);
	// System.out.println("\n=========================DictAndRLECompressed = \n" + dictAndRLECompressed);
	// System.out.println("\n=========================DictAndRLEUnCompressed = \n" + dictAndRLEUnCompressed);
	// System.out.println("\n=========================Original = \n" + sequence);
	System.out.println("Original size: " + sequence.length());
	System.out.println("DictAndRLECompressed size: " + dictAndRLECompressed.length());
	System.out.println
	    ("\nSize diff rLE to orig:   "		+ Integer.toString(sequence.length() - rLECompressed.length()));
	System.out.println
	    ("Size diff dictAndRLE to orig:   "	+ Integer.toString(sequence.length() - dictAndRLECompressed.length()));
	System.out.println(Integer.toString(sequence.length()));
	System.out.println(Integer.toString(dictAndRLECompressed.length()));

	// try(BufferedWriter writer = new BufferedWriter(new FileWriter("sequence.txt")))
	//     {
	// 	writer.write(sequence);
	//     }
	// catch (IOException e)
	//     {
	// 	e.printStackTrace();
	//     }
	// try(BufferedWriter writer = new BufferedWriter(new FileWriter("rLEAndDict.txt")))
	//     {
	// 	writer.write(dictAndRLECompressed);
	//     }
	// catch (IOException e)
	//     {
	// 	e.printStackTrace();
	//     }
	
	System.out.print("RLE: ");
	if(sequence.equals(rLEUncompressed))
	    {
		System.out.println("Equal!");
	    }
	else
	    {
		System.out.println("Not equal!");
	    }
	
	System.out.print("Dict and RLE: ");
	if(sequence.equals(dictAndRLEUnCompressed))
	    {
		System.out.println("Equal!");
	    }
	else
	    {
		for(int iter = 0; iter < sequence.length(); iter++)
		    {
			if(sequence.charAt(iter) != dictAndRLEUnCompressed.charAt(iter)) 
			    {
				System.out.println
				    ("Differ at " + iter + ", chars = " + sequence.charAt(iter) +
				     sequence.charAt(iter + 1) + sequence.charAt(iter + 2) + sequence.charAt(iter + 3) +
				     sequence.charAt(iter + 4) + sequence.charAt(iter + 5));
			    }
		    }
		System.out.println("Not equal!");
	    }
    }
}
