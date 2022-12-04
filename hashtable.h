typedef struct LinkedList LinkedList;
typedef struct Ht_item Ht_item;
typedef struct HashTable HashTable;

int hash_function(int id, int capacity);

LinkedList* allocate_list ();

LinkedList* linkedlist_insert(LinkedList* list, Ht_item* item);

void free_linkedlist(LinkedList* list); 

Ht_item* create_item(int key, int value); 

LinkedList** create_overflow_buckets(HashTable* table);

void free_overflow_buckets(HashTable* table); 

HashTable* create_table(int size); 

void free_item(Ht_item* item); 

void free_table(HashTable* table); 

void handle_collision(HashTable* table, int index, Ht_item* item);

void ht_insert(HashTable* table, int key, int value);

int ht_search(HashTable* table, int key);
