#ifndef PUNTO_D_GRAFO_H
#define PUNTO_D_GRAFO_H

#include <iostream>
#include <string>
#include <cmath>
#include <iomanip>
#include <vector>


// Estructura para los nodos de la lista enlazada (representa un lado con peso entre vértices).
struct ListaLados {
    int x;          // Vértice de origen
    int y;          // Vértice de destino
    double peso;    // Peso asociado al lado (calculado mediante distancia)
    ListaLados* next; // Puntero al siguiente nodo en la lista enlazada

    // Constructor para inicializar un nuevo lado del grafo
    ListaLados(int _x, int _y, double _peso) {
        x = _x;
        y = _y;
        peso = _peso;
        next = nullptr;
    }
};


// Clase principal del Grafo Ponderado basada en listas de adyacencia.
class Grafo {
protected:
    int numvertices;     // Número total de vértices en el grafo
    bool dirigido;       // Bandera que indica si el grafo es dirigido (true) o no dirigido (false)
    ListaLados** adjlist; // Arreglo de punteros que conforma la lista de adyacencia

public:

    // Constructor por defecto: Inicializa un grafo vacío.
    Grafo() {
        numvertices = 0;
        dirigido = false;
        adjlist = nullptr;
    }


    // Constructor con parámetros: Define el número de vértices y la direccionalidad del grafo.
    Grafo(int n, bool esdirigido) {
        numvertices = n;
        dirigido = esdirigido;
        adjlist = new ListaLados*[n];
        for (int i = 0; i < n; i++) {
            adjlist[i] = nullptr;
        }
    }


    // Destructor: Libera toda la memoria dinámica de las listas enlazadas y del arreglo de adyacencia.
    ~Grafo() {
        for (int i = 0; i < numvertices; i++) {
            ListaLados* curr = adjlist[i];
            while (curr != nullptr) {
                ListaLados* temp = curr;
                curr = curr->next;
                delete temp;
            }
        }
        delete[] adjlist;
    }

    // Método para agregar un lado con su respectivo peso calculado entre dos vértices.
    void agregarLado(int x, int y, double peso) {
        ListaLados* curr = new ListaLados(x, y, peso);
        curr->next = adjlist[x];
        adjlist[x] = curr;

        // Si el grafo es no dirigido, se registra también la arista inversa simétrica
        if (!dirigido) {
            ListaLados* currRev = new ListaLados(y, x, peso);
            currRev->next = adjlist[y];
            adjlist[y] = currRev;
        }
    }

    /**
     * Método para imprimir la estructura general del grafo consolidando los pesos
     * en formato de Matriz de Adyacencia (16x16 para 2-mers).
     */
    void imprimir() {
        // Se inicializa una matriz de 16x16 en ceros para consolidar los pesos acumulados
        std::vector<std::vector<double>> matriz(numvertices, std::vector<double>(numvertices, 0.0));

        for (int i = 0; i < numvertices; i++) {
            ListaLados* curr = adjlist[i];
            while (curr != nullptr) {
                matriz[curr->x][curr->y] += curr->peso; // Se acumula el peso en la celda correspondiente
                curr = curr->next;
            }
        }

        // Se definen las etiquetas de los 16 2-mers posibles combinando A, C, G, T
        char nucs[4] = {'A', 'C', 'G', 'T'};
        std::vector<std::string> etiquetas;
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                std::string s = "";
                s += nucs[i];
                s += nucs[j];
                etiquetas.push_back(s);
            }
        }

        // Impresión de la cabecera de la matriz con las etiquetas de columnas
        std::cout << "\n--- MATRIZ DE ADYACENCIA PONDERADA (16 x 16) ---\n\n";
        std::cout << std::setw(6) << " ";
        for (int j = 0; j < numvertices; ++j) {
            std::cout << std::setw(8) << etiquetas[j];
        }
        std::cout << "\n";

        // Impresión de las filas con sus respectivas etiquetas y valores numéricos formateados
        for (int i = 0; i < numvertices; ++i) {
            std::cout << std::setw(6) << etiquetas[i];
            for (int j = 0; j < numvertices; ++j) {
                std::cout << std::setw(8) << std::fixed << std::setprecision(4) << matriz[i][j];
            }
            std::cout << "\n";
        }
    }

    // --- MÉTODOS DE APOYO PARA GENÓMICA (2-MERS) ---


    // Convierte un carácter de nucleótido a su valor numérico correspondiente (0 a 3).

    int charToVal(char c) {
        switch (c) {
            case 'A': case 'a': return 0;
            case 'C': case 'c': return 1;
            case 'G': case 'g': return 2;
            case 'T': case 't': return 3;
            default: return -1;
        }
    }


    // Mapea un dinucleótido (2-mer) a un índice único entre 0 y 15 para estructurar los vértices del grafo.

    int getIndex2Mer(const std::string& kmer) {
        if (kmer.length() != 2) return -1;
        int v1 = charToVal(kmer[0]);
        int v2 = charToVal(kmer[1]);
        if (v1 == -1 || v2 == -1) return -1;
        return (v1 * 4) + v2;
    }
};

#endif
