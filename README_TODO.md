# Strategy for the new mapper:

The new mapper should store the ids of a fragments vector linking to the sequences.
The sequences on the other hand are of the trait rustody::traits::BinaryMatcher and implement as data structure a Vec\<u8\>. The new fast_mapper.map entries need therefore be a (\<genes_id: usize\> and \<start: usize\> ). The genes vector will be filled with the new BiaryMatcher objects and we then can simply match the whole read against this gene.

## first objective:

Create the GenesData class and implement BinaryMatcher for that. I also need IntToStr's seq_at_position function, but with an start and stop value. So a name like subset or slice would make more sens.
