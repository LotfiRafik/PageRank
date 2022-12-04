#include "stdio.h"
#include "stdlib.h"
#include "hashtable.h"

// Define the Hash Table Item here
struct Ht_item {
    int key;
    int value;
};

// Define the Hash Table here
struct HashTable {
    // Contains an array of pointers
    // to items
    Ht_item** items;
    LinkedList** overflow_buckets;
    int size;
    int count;
};

// Define the Linkedlist here
struct LinkedList {
    Ht_item* item; 
    LinkedList* next;
};

int hash_function(int id, int capacity) {
    return id % capacity;
}

LinkedList* allocate_list () {
    // Allocates memory for a Linkedlist pointer
    LinkedList* list = (LinkedList*) malloc (sizeof(LinkedList));
    return list;
}

LinkedList* linkedlist_insert(LinkedList* list, Ht_item* item) {
    // Inserts the item onto the Linked List
    if (!list) {
        LinkedList* head = allocate_list();
        head->item = item;
        head->next = NULL;
        list = head;
        return list;
    } 
    
    else if (list->next == NULL) {
        LinkedList* node = allocate_list();
        node->item = item;
        node->next = NULL;
        list->next = node;
        return list;
    }

    LinkedList* temp = list;
    while (temp->next->next) {
        temp = temp->next;
    }
    
    LinkedList* node = allocate_list();
    node->item = item;
    node->next = NULL;
    temp->next = node;
    
    return list;
}

void free_linkedlist(LinkedList* list) {
    LinkedList* temp = list;
    while (list) {
        temp = list;
        list = list->next;
        free(temp->item);
        free(temp);
    }
}

Ht_item* create_item(int key, int value) {
    // Creates a pointer to a new hash table item
    Ht_item* item = (Ht_item*) malloc (sizeof(Ht_item));
    item->key = key;
    item->value = value;

    return item;
}

LinkedList** create_overflow_buckets(HashTable* table) {
    // Create the overflow buckets; an array of linkedlists
    LinkedList** buckets = (LinkedList**) calloc (table->size, sizeof(LinkedList*));
    for (int i=0; i<table->size; i++)
                   buckets[i] = NULL;
    return buckets;
}

void free_overflow_buckets(HashTable* table) {
    // Free all the overflow bucket lists
    LinkedList** buckets = table->overflow_buckets;
    for (int i=0; i<table->size; i++)
          free_linkedlist(buckets[i]);
    free(buckets);
}

HashTable* create_table(int size) {
    // Creates a new HashTable
    HashTable* table = (HashTable*) malloc (sizeof(HashTable));
    table->size = size;
    table->count = 0;
    table->items = (Ht_item**) calloc (table->size, sizeof(Ht_item*));
    for (int i=0; i<table->size; i++)
        table->items[i] = NULL;
    table->overflow_buckets = create_overflow_buckets(table);

    return table;
}

void free_item(Ht_item* item) {
    // Frees an item
    free(item);
}

void free_table(HashTable* table) {
    // Frees the table
    for (int i=0; i<table->size; i++) {
        Ht_item* item = table->items[i];
        if (item != NULL)
            free_item(item);
    }

    // Free the overflow bucket linked linkedlist and it's items
    free_overflow_buckets(table);
    free(table->items);
    free(table);
}

void handle_collision(HashTable* table, int index, Ht_item* item) {
    LinkedList* head = table->overflow_buckets[index];

    if (head == NULL) {
        // We need to create the list
        head = allocate_list();
        head->item = item;
        table->overflow_buckets[index] = head;
        return;
    }
    else {
        // Insert to the list
        table->overflow_buckets[index] = linkedlist_insert(head, item);
        return;
    }
 }

void ht_insert(HashTable* table, int key, int value) {
	// Create the item
	Ht_item* item = create_item(key, value);

    // Compute the index    
    int index = hash_function(key, table->size);

	Ht_item* current_item = table->items[index];
	
	if (current_item == NULL) {
	    // Key does not exist.
	    if (table->count == table->size) {
	        // Hash Table Full
	        printf("Insert Error: Hash Table is full\n");
	        return;
	    }
	    
	    // Insert directly
	    table->items[index] = item; 
	    table->count++;
	}

	else {
    		// Scenario 1: We only need to update value
    		if (current_item->key == key) {
        		table->items[index]->value = value;
        		return;
    		}
    
    	else {
			// Scenario 2: Collision
		  	// We will handle case this a bit later
          	handle_collision(table, index, item);
          	return;
    	}
	}
}


int ht_search(HashTable* table, int key) {
    // Searches the key in the hashtable
    // and returns NULL if it doesn't exist
    int index = hash_function(key, table->size);
    Ht_item* item = table->items[index];
    LinkedList* head = table->overflow_buckets[index];

    // Ensure that we move to items which are not NULL
    while (item != NULL) {
        if (item->key == key)
            return item->value;
        if (head == NULL)
            return -1;
        item = head->item;
        head = head->next;
    }
    return -1;
}

