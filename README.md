# REL And Dict Compression
## A Program To Compress Text That Uses Both RLE And Dictionary Compression On The Input Text
*this is a program that I wrote while working on another project. That project is written in C++ and I had a basic RLE implementation where it would be able to compress repeated sequences of length 1 (that is 1 character repeated n times, where n is above some base number.) I wanted to create a better compression scheme for that program. But I knew that it would be at least a little complicated, so I decided to write a separate program to flesh it out. I decided to write this program in Java since I thought this would probably speed up the development process a little and I was (and am) pretty rusty with Java, so I also thought it would be a good chance to brush up on it a little.
The program first compresses the input text using an RLE scheme that is implemented using a dynamic programming approach and can represent repeated sequences of n characters (that is to say that it can represent repeated sequences, where the sequences repeated are more than 1 character long.)*

*Then the program performs dictionary compression on the RLE encoded string. The member function that does this is passed a number which is the size of string to try to find multiple occurrences of. It uses a hash table that keeps track of the number of occurrences of a given pre-image (source string that produced the hash). This way it can more efficiently check if it has already seen a given string. Then for all strings that occur more than once it goes over all of the locations that the pre-images were found at for a given pre-image and makes sure that there are no overlaps from left to right. If there is an over lap that pre-image location is discarded. After this it makes sure that that cumulative size of the remaining pre-image locations for the given pre-image are less than the size of replacing them with references to one copy of the pre-image. If there would be a size benefit from replacing the pre-image locations with references to one copy of the pre-image, it is done. The function that does this is called multiple times with decreasing sizes for the pre-image search. It is also passed special characters to use to denote the start of a reference\/pre-image and the end of what we will call the \"dictionary area\". Then to undo the dictionary compression an un-compression function is called and passed the special characters in the opposite order to the compression function.*

*It should be noted that without modification this program isn't that useful for text compression because*
1. *It uses ASCII numbers for the length for the RLE encoding. So the text shouldn't contain numbers. This could be alleviated by making the repeat count a fixed length (but this would also make it less efficient.*
2.  *It uses Hiragana and Katakana numbers as special characters for both the RLE and dictionary encoding. This is not just a problem because the program cannot be used on Japanese text. It is a problem because when the compressed representation is written out (using a BufferedWriter in Java), it will be written out in UTF-16, which is a variable length encoding. Thus the compressed text can be less characters than the original text, but take up more space. However as far as we know it should take up less memory when the program is running because Java will always use at least two bytes for a given character (even if the UTF-16 representation is less, it will be padded.) This could be fixed by just using character above 128, but less than 256. However we are using that range for the pre-image reference numbers. Anyway all of that is to say that the program could be modified so that there would be at least some gain when written to disc.*

*It is of note that the above described problems won't be a problem for the
intended use of the compression scheme since it is intended to be used in a
situation where we have a number of 16 bit values that are always quite a bit
less than the full range of a 16 bit number. This means that we can use values
at the top of the range for special values and we can just use a fixed length
encoding for length for RLE and for pre-image references. So basically because
of the nature of the data we want to encode we should benefit from this encoding
scheme. We should also note that the RLE compression part of the program is
fairly inefficient in terms of time complexity. However that data we are going
to compress is always of a fixed (and manageable length). Also note that the
program definitely needs a bit of work, but as mentioned this is really for
another project, so we will be translating it into C++ and cleaning it up a bit
and doing much more testing.*


*We've found a bug in the code around 646 that results in less optimal compression in terms of the dictionary side of things. However we haven't fixed it, because as stated above we wrote this code to use in a C++ project and we've fixed the code there. We are just making a note here so that if we ever want to use this code for something, we know to fix it.*
