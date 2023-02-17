typedef struct LinkedList LinkedList;
typedef struct Ht_item Ht_item;
typedef struct HashTable HashTable;

// Define the Hash Table Item here
typedef struct Ht_item {
    long key;
    int value;
} Ht_item;

// Define the Linkedlist here
typedef struct LinkedList {
    Ht_item* item; 
    LinkedList* next;
} LinkedList;

// Define the Hash Table here
typedef struct HashTable {
    // Contains an array of pointers
    // to items
    Ht_item** items;
    LinkedList** overflow_buckets;
    long size;
    long count;
} HashTable;



long hash_function(long id, long capacity);

LinkedList* allocate_list ();

LinkedList* linkedlist_insert(LinkedList* list, Ht_item* item);

void free_linkedlist(LinkedList* list); 

Ht_item* create_item(long key, int value); 

LinkedList** create_overflow_buckets(HashTable* table);

void free_overflow_buckets(HashTable* table); 

HashTable* create_table(long size); 

void free_item(Ht_item* item); 

void free_table(HashTable* table); 

void handle_collision(HashTable* table, long index, Ht_item* item);

void ht_insert(HashTable* table, long key, int value);

int ht_search(HashTable* table, long key);
