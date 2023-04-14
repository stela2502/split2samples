# I am wondering how to implement a fast mapper.

And I had a discussion with ChatGPT to point me to the right solution.

The last idea I (we) had was: The initial 8bp lookup vector could be used as a quick filter to quickly narrow down the search space before doing more detailed matching using the radix_trie. The radix_trie can be optimized for memory usage by breaking it down into smaller sub-trees if necessary. By doing this, you can minimize the amount of memory used while still maintaining fast lookup times.

I need to adjust the GeneIDs class to use this idea in the back.


# The whole discussion:


Me:
Can you tell me how telephone numbers are stored in Rust?

In Rust, telephone numbers can be stored using different data types depending on the requirements of the application. One common way to represent a telephone number is by using a string data type. This can be done using Rust's built-in String type, which provides dynamic string storage and manipulation.

For example, a telephone number can be stored in Rust like this:

rust

let telephone_number = String::from("123-456-7890");

This stores the telephone number "123-456-7890" in a String variable called telephone_number.

Alternatively,

Me:
OK that was too simple a question. I want to store ALL telephone numbers and the aim of the tool should be to identify information stored with the numbers lightning fast.

Storing all possible telephone numbers in memory would not be a feasible approach due to the large number of possible phone numbers. Instead, you would need to use some form of indexing or hashing to efficiently search and retrieve information associated with phone numbers.

One commonly used approach is to use a database management system (DBMS) such as SQLite, MySQL or PostgreSQL to store the phone numbers and associated information. These systems provide efficient indexing and querying capabilities, allowing lightning fast retrieval of information associated with a phone number.

When storing phone numbers in a database, you can choose an appropriate data type to represent the phone number. One common approach is to store the phone number as a string, similar to the example I provided earlier. However, you can also store the phone number as a numeric value using Rust's built-in integer types such as u64 or i64. This can be more efficient for some operations, but may require additional handling to format the phone number for display purposes.

To optimize searching and retrieval, you can create an index on the phone number column in the database. This will allow the DBMS to efficiently locate the records associated with a particular phone number, even in a large dataset. Additionally, you can consider using caching mechanisms or other optimizations to further improve performance.
Me:
I once heard that a developer significantly improved mapping times if sequences were stored binary like telephone numbers. I am thinking about a u32 stored in a Vec<u32>, but use the 16 bp sequence as index in this Vec. How large would a Vec<u32> be if it would be the maximum length?

If you are storing a binary sequence of 16 base pairs as an index in a Vec<u32>, and assuming you are using one base pair per two bits (since there are four possible values for each base pair, and two bits can represent four possible values), you would need a total of 32 bits (or 4 bytes) to represent each index.

The maximum number of 16-base pair sequences is 4^16, or 65,536. If you were to store each of these as an index in a Vec<u32>, you would need a vector of size 65,536 * 4 bytes, which is equal to 262,144 bytes or 256 kilobytes.

Note that this calculation only considers the space required to store the indices themselves, and not the space required to store the associated data or any overhead associated with the vector itself. Also, if you plan to use this data structure to map sequences to phone numbers, you would need to ensure that each sequence maps to a unique phone number or set of phone numbers to avoid collisions.
Me:
Oh - that is actually not a lot of memory. But this gets singificantly bigger if I use u64 - or?

Yes, if you use u64 instead of u32, the size of your vector will increase. Each u64 value requires 8 bytes of memory, so if you use u64 to store the indices, the size of your vector would be:

python

65,536 * 8 bytes = 524,288 bytes = 512 kilobytes

As you can see, the size of the vector increases significantly when using u64 instead of u32. However, 512 kilobytes is still a relatively small amount of memory by modern standards, and may be acceptable depending on the requirements of your application.

Note that if you are planning to store additional data associated with each index (such as phone numbers), the overall memory requirements will be higher. In this case, you would need to consider the space required for the index as well as the space required for the associated data.
Me:
That is true. So if I now use the index as u64, too it would become even bigger.

Yes, if you use u64 for both the index and the phone number, the overall memory usage will increase further.

Assuming you are using u64 for both the index and the phone number, the size of your data structure for storing all possible phone numbers would be:

python

65,536 * (8 bytes + 8 bytes) = 1,048,576 bytes = 1 megabyte

As you can see, the size of the data structure increases to 1 megabyte when using u64 for both the index and the phone number. This is still a relatively small amount of memory by modern standards, but it may be significant depending on the requirements of your application.

Note that this calculation assumes that you are storing only a single phone number for each index. If you are storing multiple phone numbers for each index, the size of your data structure will be larger. Additionally, this calculation does not include any overhead associated with the data structure itself or any additional data that you may need to store.
Me:
I am not buying your answer ;-). Rust tells me the max value of a u64 is 18_446_744_073_709_551_615 so your initial 65,536 is incorrect - or?

You are correct, I apologize for the mistake in my earlier response. The maximum value of a u64 is indeed 18,446,744,073,709,551,615, which is much larger than the value I used in my previous calculations.

Assuming you are using a u64 for both the index and the phone number, the size of your data structure for storing all possible phone numbers would be:

python

4^16 * (8 bytes + 8 bytes) = 4,398,046,511,104 bytes = 4 terabytes

As you can see, the size of the data structure increases to 4 terabytes when using u64 for both the index and the phone number. This is a significant amount of memory and may be impractical for most systems.

However, note that storing all possible phone numbers is not a realistic use case, as most phone numbers are not in use and do not need to be stored. Instead, you would only need to store phone numbers that are actually in use, which would significantly reduce the memory requirements.
Me:
Right - so now I see why my idea is doomed. So if I scale that back to u32 - please run this calculation instead...

Sure! If you are using u32 for both the index and the phone number, the size of your data structure for storing all possible phone numbers would be:

python

4^16 * (4 bytes + 4 bytes) = 1,099,511,627,776 bytes = 1 terabyte

As you can see, the size of the data structure is still quite large at 1 terabyte, but it is more practical than using u64. However, as I mentioned earlier, storing all possible phone numbers is not a realistic use case, so the actual memory requirements would be much lower.
Me:
Thought you would fall into that trap. Your 4^16 is wrong.

You're right, I apologize for the error. The correct value for the number of possible 16-bit sequences is 65,536, which is equal to 2^16.

So if you are using u32 for both the index and the phone number, the size of your data structure for storing all possible phone numbers would be:

python

65,536 * (4 bytes + 4 bytes) = 524,288 bytes = 512 kilobytes

As you can see, the size of the data structure is much smaller when using u32 instead of u64. This size may be acceptable depending on the requirements of your application. However, keep in mind that this calculation assumes that you are only storing a single phone number for each index, and the actual memory requirements may be higher if you need to store additional data with each index.
Me:
I feel like I need to tell you more than you help me :-D: U32 MAX value is 4_294_967_295. 

You are correct, I apologize for the mistake in my earlier response. The maximum value of a u32 is indeed 4,294,967,295.

Assuming you are using u32 for both the index and the phone number, the size of your data structure for storing all possible phone numbers would be:

python

65,536 * (4 bytes + 4 bytes) = 524,288 bytes = 512 kilobytes

As you can see, the size of the data structure is much smaller than when using u64. However, keep in mind that this calculation assumes that you are only storing a single phone number for each index, and the actual memory requirements may be higher if you need to store additional data with each index.
Me:
 4,294,967,295 * (32 bytes) =

Multiplying 4,294,967,295 by 32 bytes gives:

python

4,294,967,295 * 32 bytes = 137,438,953,440 bytes = 128 gigabytes

As you can see, using a u32 to store each phone number would result in a significantly larger data structure if you need to store all possible phone numbers. This is why it's important to consider the actual requirements of your application and only store the phone numbers that you need to work with.
Me:
But storing them all would significantly increase the speed of the application. Its always a question if speed should be prioritized. So if I would boil it down further and only look for u8 (4bp) I would be at the 2^8 * 64 bytes if I use pointers to the possible next 4bp parts. Stacking this logics up to 4 iterations would end me with 16 bp again. Assuming I have possibly only 10% of the initial 8 byte vector filled and each other has only like 10 entries - how much memory would I need?

If I understand correctly, you are proposing a trie-like data structure where each node has 256 pointers to the possible next 4-bit parts, and you are only interested in 16-bit sequences. In this case, the size of the data structure would be:

python

256^(16/4) * 8 bytes per pointer = 1,208,925,819,614,629,174,706,176 bytes = 1.2 x 10^24 bytes

This is an incredibly large amount of memory, much larger than what would be practical to use on a modern computer. Keep in mind that this calculation assumes that every possible 16-bit sequence is stored in the trie, and even if only 10% of the initial 8-byte vector is filled, the resulting data structure would still be many orders of magnitude larger than what is practical.

It's important to consider the actual requirements of your application and balance the tradeoffs between speed and memory usage. Storing all possible phone numbers is likely not practical, but storing a subset of phone numbers that are frequently accessed could be a good compromise.
Me:
What is the log2 of 65,536

The logarithm base 2 of 65,536 is 16, because 2^16 = 65,536.

So mathematically, we can write:

log2(65,536) = 16
Me:
so 16bit index and u64 entry would work.

Yes, using a 16-bit index and a 64-bit entry would work for storing telephone numbers. The data structure would be able to store all possible 16-bit sequences as indices, and a 64-bit entry for each index would allow for storing additional information associated with each phone number.

Assuming that you have additional data to store with each phone number, the size of the data structure would be:

python

65,536 * (8 bytes + additional data size) bytes

The additional data size would depend on what other information you need to store with each phone number. However, keep in mind that even with a 64-bit entry, storing all possible phone numbers may not be practical depending on the requirements of your application.
Me:
u16 and the entry is a pointer to another object - but these are most of the time empty. How much would just the vector of pointers be?
2 / 2

If you're using a Vec of pointers to another object to store the telephone numbers, the size of the vector would depend on how many telephone numbers you need to store.

Assuming you need to store all possible 16-bit sequences, the size of the vector would be:

arduino

65,536 * size_of::<*const T>()

where T is the type of the object that the pointers are pointing to. size_of::<*const T>() gives the size of a pointer to type T in bytes.

For example, if you're storing pointers to objects of type MyObject, and assuming that size_of::<*const MyObject>() is 8 bytes on your system, the size of the vector would be:

python

65,536 * 8 bytes = 524,288 bytes = 512 kilobytes

Note that this calculation only takes into account the size of the vector itself and not the size of the objects being pointed to. If the objects being pointed to are large, this could add significant memory usage to your application.
Me:
Is an 8 bytes pointer common in Rust=

Yes, an 8-byte pointer is common in Rust for 64-bit systems. This is because pointers are typically implemented as 64-bit memory addresses on these systems.

However, the exact size of a pointer can depend on the specific system and the implementation of Rust that you're using. For example, on a 32-bit system, pointers would typically be 4 bytes in size. Additionally, some Rust implementations or platforms may use different pointer sizes for specific purposes or in certain contexts.
Me:
That could be feasible as I use this 512 kb object just for the really fast initial lookup. So if I would not use 8 bp initial match but 16 bp we are 4,294,967,295 * 8 bytes. How much is that?

If you're storing all possible 16-bit sequences as keys in a HashMap or similar data structure, then the size of the data structure would depend on the size of the values associated with each key.

Assuming that you're using 8-byte values (e.g., pointers) to associate with each key, then the size of the data structure would be:

python

4,294,967,295 * 8 bytes = 34,359,738,360 bytes = 32 gigabytes

This is a very large amount of memory, and it's unlikely that you would be able to store all possible phone numbers in memory in this way. Additionally, keep in mind that using a HashMap or similar data structure can have additional overhead beyond just the memory usage, such as the cost of hashing and collision resolution.
Me:
Now I need to really blow your calculations: If I use 100bp of every human genes - let's say 40000 times 100 bp. Repeat sequences are not included and I use this strategy to use three possible 8bp sequences per gene - this should reduce the number of sequences for the first 40,000 significantly for the second step. So I would not need a memory inensive second step - right? 

It's true that if you can reduce the number of possible sequences to be stored in memory, you can reduce the memory requirements for your program.

However, keep in mind that if you're using a hashing-based data structure like a HashMap, the time required to search for a key may depend on the number of keys in the data structure, even if many of them are empty. This is because the hash function still needs to be applied to each key to determine its index in the data structure, and the index needs to be searched for the key even if the value associated with the key is empty.

Therefore, even if you reduce the number of possible sequences for the first 40,000 significantly, you may still need to use a data structure that can handle a large number of keys efficiently, especially if you need to perform fast searches across all of the possible keys.

One possible solution could be to use a Trie data structure to store the sequences, which can be more memory-efficient than a HashMap for certain types of data. A Trie can also be searched efficiently for all keys that match a particular prefix, which could be useful for your use case.
Me:
What is a Trie data structure?

A Trie (pronounced "try") data structure is a tree-like data structure used for efficient retrieval of keys that match a given prefix. It is also known as a digital tree, radix tree, or prefix tree.

In a Trie, each node represents a prefix of one or more keys. The root node represents an empty prefix. Each edge leaving a node represents a character, and the concatenation of the characters along a path from the root to a leaf node represents a complete key in the Trie.

To insert a key into a Trie, the key is broken down into its individual characters, and a path is created in the Trie for each character in the key. The last node in the path is marked as a leaf node, and the value associated with the key is stored at that node.

To search for a key in a Trie, the search key is broken down into its individual characters, and the Trie is traversed one character at a time. If a path corresponding to the search key exists in the Trie, then the value associated with the key can be found at the leaf node at the end of the path.

One of the main advantages of a Trie is that it can perform prefix searches efficiently. A prefix search involves finding all keys in the Trie that start with a given prefix. This can be done by traversing the Trie to the node representing the prefix, and then traversing all the subtree below that node to find all leaf nodes.

Tries can be implemented using arrays, linked lists, or hash tables, among other data structures. In Rust, the trie and radix_trie crates provide implementations of Tries that you can use.
Me:
Possibly the thing I want. What is the difference between trie and radix_trie

Radix trie is a variant of the trie data structure where nodes that have only one child node are merged with their child. This allows for space optimization since fewer nodes are needed to represent the same set of keys. For example, if a trie has the keys "car" and "card", the trie would have four nodes: the root, "c", "a", and "r", with "d" and "card" as children of "r". In a radix trie, the "r" and "d" nodes would be merged into a single node, reducing the number of nodes to three.

The radix trie is a space-optimized variant of the trie data structure, and is particularly useful in cases where the keys are long and have common prefixes. However, it may have slightly slower lookup times than a regular trie due to the need for additional processing to merge nodes.

In summary, the main difference between a trie and a radix trie is the way they handle nodes with a single child.
Me:
Is this a recursive approach? Say I have u64 I want to store and it would make sense to create a binary radix_trie from that as the probability after 32 bits would almost be zero that there is more than one leaf. Or would it be more efficient to use a two bit 'word' for the trie?

Yes, both trie and radix trie data structures use a recursive approach. In a trie, each node represents a character or a part of a key and has pointers to child nodes representing subsequent characters or parts of the key. In a radix trie, nodes are merged if they have only one child, so the structure becomes more compact.

Regarding your second question, if you are dealing with 64-bit values, a binary radix trie might not be the best option as it would require a large number of levels to represent all possible values. Using a two-bit word would be more efficient in terms of space, but it would also increase the height of the tree and possibly slow down lookups.

In general, the best data structure to use depends on the specific requirements of your application, such as the size of the data set, the frequency of insertions and lookups, and the memory constraints of your system.
Me:
Thank you. Say I have 120Gb of memory and the frequency of insertions is zero as after creation it would be a look up tree.

If you have 120GB of memory available and the frequency of insertions is zero, you can use a radix trie to store your data. With a radix trie, you can store your data in a compact way that minimizes memory usage and provides fast lookups.

However, keep in mind that the memory usage of a radix trie depends on the specific characteristics of your data set. For example, if your data is evenly distributed and there are no long common prefixes, a radix trie might be more memory efficient than a standard trie. But if your data has long common prefixes, a radix trie might not be the best option, as it can result in a large number of nodes with only one child.

In any case, with 120GB of memory, you should be able to store a large amount of data in a radix trie and still have room for other data structures and system overhead.
Me:
So to bring this to a max. If I use the initial 8bp lookup vector with 512 kilobytes of size and store a set of pointers of radix_trie objects in this vector I would have rather small radix_trie objects and a really fast match of the initial 8bp sequences.

Yes, that sounds like a reasonable approach. The initial 8bp lookup vector could be used as a quick filter to quickly narrow down the search space before doing more detailed matching using the radix_trie. The radix_trie can be optimized for memory usage by breaking it down into smaller sub-trees if necessary. By doing this, you can minimize the amount of memory used while still maintaining fast lookup times.
Me:
Nice - I'll try that.

Great, I hope it works well for your use case! Don't hesitate to reach out if you have any further questions.