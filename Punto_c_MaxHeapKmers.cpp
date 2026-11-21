#ifndef PUNTO_C_MAX_HEAP_KMERS_H
#define PUNTO_C_MAX_HEAP_KMERS_H

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>

// Nota: KmerFrequency ya está definido en Punto_b_doubleHashing.cpp,
// no se vuelve a declarar para que no ocurran errores de definición.


// Clase principal que gestiona la estructura de Max-Heap orientada a frecuencias de K-mers.
class MaxHeapKmers {
private:
    std::vector<KmerFrequency> heap; // Vector interno que representa el árbol binario del montículo

    // Se usa para ordenar el heap despues de una insercion (percolate-up)
    void heapifyUp(int index) {
        while (index > 0) {
            int parent = (index - 1) / 2;
            // Ordenamos de mayor a menor frecuencia (el mayor valor asciende a la raíz)
            // Mientras el padre sea MENOR que el actual, se intercambian
            if (heap[index].frequency > heap[parent].frequency) {
                std::swap(heap[index], heap[parent]);
                index = parent;
            } else {
                break;
            }
        }
    }

    // Se usa para ordenar el heap despues de una eliminacion (percolate-down)
    void heapifyDown(int index) {
        int size = heap.size();
        while (2 * index + 1 < size) {
            int left = 2 * index + 1;
            int right = 2 * index + 2;
            int largest = index;

            // Comparamos con el hijo izquierdo
            if (left < size && heap[left].frequency > heap[largest].frequency) {
                largest = left;
            }
            // Comparamos con el hijo derecho
            if (right < size && heap[right].frequency > heap[largest].frequency) {
                largest = right;
            }

            // Si el mayor no es el nodo actual, intercambiamos y continuamos bajando
            if (largest != index) {
                std::swap(heap[index], heap[largest]);
                index = largest;
            } else {
                break;
            }
        }
    }

public:
    //Constructor por defecto: Inicializa un Max-Heap vacío.
    MaxHeapKmers() {}

    /**
     * Construye el Max-Heap de manera óptima (O(N)) a partir de un vector
     * de K-mers y frecuencias obtenidos de la tabla hash.
     */
    void buildHeap(const std::vector<KmerFrequency>& datos) {
        heap = datos;
        // Se aplica heapifyDown desde el último nodo no hoja hasta la raíz
        for (int i = (heap.size() / 2) - 1; i >= 0; --i) {
            heapifyDown(i);
        }
    }

    /**
     * Extrae y retorna el elemento con mayor frecuencia (la raíz del Max-Heap),
     * reorganizando el montículo subsiguiente.
     */
    KmerFrequency extractMax() {
        if (heap.empty()) {
            return {"", -1};
        }
        KmerFrequency maxElement = heap[0]; // La raíz contiene el valor máximo
        heap[0] = heap.back();
        heap.pop_back();
        if (!heap.empty()) {
            heapifyDown(0);
        }
        return maxElement;
    }

    /**
     * Imprime en consola los 'N' K-mers más frecuentes utilizando una copia temporal
     * del heap para no alterar la estructura original durante las extracciones.
     */
    void printTopN(int n) {
        std::cout << "\n--- TOP " << n << " K-MERS MAS FRECUENTES (MAX-HEAP) ---\n";

        // Creamos una copia del heap para preservar los datos originales intactos
        std::vector<KmerFrequency> copiaHeap = heap;

        // Reconstruimos el heap localmente con la copia
        MaxHeapKmers tempHeap;
        tempHeap.buildHeap(copiaHeap);

        int count = 0;
        while (count < n && !tempHeap.heap.empty()) {
            KmerFrequency kf = tempHeap.extractMax();
            std::cout << count + 1 << ". K-mer: " << kf.kmer
                      << " | Frecuencia: " << kf.frequency << "\n";
            count++;
        }
    }
};

#endif
