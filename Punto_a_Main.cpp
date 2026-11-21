#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>

// Inclusión de módulos desarrollados para las diferentes estructuras de datos
#include "Punto_b_hashLinearProbing.cpp"
#include "Punto_b_doubleHashing.cpp"
#include "Punto_c_MaxHeapKmers.cpp"
#include "Punto_d_Grafo.cpp"

using namespace std;

/**
 * Estructura para representar una secuencia de ADN individual.
 * Almacena el identificador único (ID) y la cadena de nucleótidos asociada.
 */
struct DNASequence {
    string id;
    string sequence;
};

/**
 * Función encargada de leer, parsear y sanitizar archivos en formato FASTA.
 * - Abre el archivo de texto especificado.
 * - Filtra saltos de línea adicionales y líneas vacías.
 * - Extrae el ID de la cabecera (marcada con '>') y concatena los bloques de secuencia.
 * - Estandariza los nucleótidos a mayúsculas para asegurar consistencia.
 */
void cargarArchivoFasta(vector<DNASequence*>& listaSecuencias, const string& nombreArchivo) {
    // Abre el archivo de FASTA especificado
    ifstream archivo(nombreArchivo);

    // Validación de apertura del archivo
    if (!archivo.is_open()) {
        cout << "[ERROR] No se pudo abrir el archivo FASTA: " << nombreArchivo << "\n";
        return;
    }

    string linea;
    string idActual = "";
    string secuenciaActual = "";

    // Lectura línea por línea del archivo FASTA
    while (getline(archivo, linea)) {
        // Manejo de compatibilidad de saltos de carro en diferentes sistemas operativos (\r)
        if (!linea.empty() && linea.back() == '\r') {
            linea.pop_back();
        }

        // Ignorar líneas vacías
        if (linea.empty()) continue;

        // Detección de una nueva cabecera FASTA (comienza con '>')
        if (linea[0] == '>') {
            if (!idActual.empty()) {
                // Almacenamiento de la secuencia previa antes de cambiar de cabecera
                DNASequence* nuevaSeq = new DNASequence();
                nuevaSeq->id = idActual;
                nuevaSeq->sequence = secuenciaActual;
                listaSecuencias.push_back(nuevaSeq);
                secuenciaActual = "";
            }
            // Extracción del identificador único de la línea de cabecera
            stringstream ss(linea.substr(1));
            ss >> idActual;
        } else {
            // Acumulación de los fragmentos de la secuencia de ADN convirtiéndolos a mayúsculas
            for (char &c : linea) {
                c = toupper(c);
            }
            secuenciaActual += linea;
        }
    }

    // Registro de la última secuencia procesada al finalizar el archivo
    if (!idActual.empty()) {
        DNASequence* ultimaSeq = new DNASequence();
        ultimaSeq->id = idActual;
        ultimaSeq->sequence = secuenciaActual;
        listaSecuencias.push_back(ultimaSeq);
    }

    archivo.close();
    cout << "Se han cargado exitosamente " << listaSecuencias.size() << " secuencias desde el archivo FASTA.\n";
}

/**
 * Función para extraer subcadenas de longitud 'n' (K-mers) de una secuencia de ADN
 * y alimentar la Tabla Hash basada en Double Hashing para conteo de frecuencias.
 */
void extraerYKmers(const string& secuenciaDNA, int n, DoubleHashTable& tablaKmers) {
    if (secuenciaDNA.length() < (size_t)n) return;

    // Recorrido de la cadena para extraer cada k-mer deslizante
    for (size_t i = 0; i <= secuenciaDNA.length() - n; ++i) {
        string kmer = secuenciaDNA.substr(i, n);
        tablaKmers.insertOrIncrement(kmer);
    }
}

int main() {
    // Vector principal en memoria dinámica para almacenar las secuencias de ADN cargadas
    vector<DNASequence*> listaSecuencias;

    // 1. Carga automática inicial del archivo FASTA al arrancar el programa
    cout << "Cargando archivo FASTA...\n";
    cargarArchivoFasta(listaSecuencias, "HIV2_ALL_2021_env_DNA.2000.fasta");

    // Inicialización y llenado automático de la Tabla Hash 1 (Linear Probing) con los IDs y secuencias
    HashTable tablaSecuencias1(10);
    for (DNASequence* seq : listaSecuencias) {
        tablaSecuencias1.insert(seq->id, seq->sequence);
    }

    int opcion;
    string idBusqueda;
    int nKmer = 3; // Longitud por defecto para los K-mers (ej. Trímeros o 3-mers)

    // Bucle principal del menú interactivo de la aplicación
    do {
        cout << "\n===================================================\n";
        cout << "      MINI SISTEMA DE ANALISIS GENOMICO (HIV)      \n";
        cout << "===================================================\n";
        cout << "1. Ver resumen de secuencias de ADN cargadas\n";
        cout << "2. Tabla Hash 1: Buscar secuencia por ID (Linear Probing)\n";
        cout << "3. Tabla Hash 2: Extraer y contar K-mers (Double Hashing)\n";
        cout << "4. Max-Heap: Ver patrones (K-mers) mas frecuentes\n";
        cout << "5. Grafo: Generar matriz de adyacencia (2-mers)\n";
        cout << "6. Salir\n";
        cout << "Seleccione una opcion: ";

        // Validación de entrada numérica en el menú
        if (!(cin >> opcion)) {
            cout << "\n[ERROR] Entrada invalida. Por favor, ingrese un NUMERO del 1 al 6.\n";
            cin.clear();
            cin.ignore(10000, '\n');
            continue;
        }
        cin.ignore();

        switch(opcion) {
            case 1: { // Opción 1: Listar resumen de IDs cargados
                cout << "\n--- LISTA DE SECUENCIAS CARGADAS ---\n";
                int limite = (listaSecuencias.size() < 220) ? listaSecuencias.size() : 220;
                for (int i = 0; i < limite; ++i) {
                    cout << "ID: " << listaSecuencias[i]->id << "\n";
                }
                cout << "Total de secuencias en memoria: " << listaSecuencias.size() << "\n";
                break;
            }
            case 2: { // Opción 2: Búsqueda rápida por ID usando Tabla Hash 1 (Linear Probing)
                cout << "\n--- BUSQUEDA EN TABLA HASH 1 (Linear Probing) ---\n";
                cout << "Ingrese el ID de la secuencia a buscar: ";
                getline(cin, idBusqueda);

                string* resultado = tablaSecuencias1.lookup(idBusqueda);
                if (resultado != nullptr) {
                    cout << "[ENCONTRADO] Secuencia asociada al ID encontrada con exito.\n";
                    cout << "-> Longitud de la secuencia: " << resultado->length() << " pares de bases\n";
                    cout << "-> Fragmento inicial: " << resultado->substr(0, 50) << "...\n";
                } else {
                    cout << "[NO ENCONTRADO] No existe una secuencia con ese identificador.\n";
                }
                break;
            }
            case 3: { // Opción 3: Extracción y conteo de K-mers con Tabla Hash 2 (Double Hashing)
                cout << "\n--- TABLA HASH 2 (Double Hashing) ---\n";
                cout << "Ingrese la longitud del K-mer (n): ";
                cin >> nKmer;
                cin.ignore();

                DoubleHashTable tablaKmers(nKmer);
                // Procesamiento de la primera secuencia como muestra representativa
                if (!listaSecuencias.empty()) {
                    extraerYKmers(listaSecuencias[0]->sequence, nKmer, tablaKmers);
                    cout << "Se han extraido y contado los " << nKmer << "-mers de la primera secuencia.\n";
                    tablaKmers.display();
                }
                break;
            }
            case 4: { // Opción 4: Visualización de patrones frecuentes ordenados mediante Max-Heap
                cout << "\n--- MAX-HEAP DE K-MERS FRECUENTES ---\n";
                if (listaSecuencias.empty()) {
                    cout << "[ERROR] No hay secuencias cargadas.\n";
                    break;
                }

                cout << "Procesando k-mers (" << nKmer << "-mers) de la secuencia actual...\n";
                DoubleHashTable tablaKmers(nKmer);
                extraerYKmers(listaSecuencias[0]->sequence, nKmer, tablaKmers);

                // Extracción de los k-mers válidos generados en la tabla hash
                vector<KmerFrequency> datosKmers = tablaKmers.obtenerKmersValidos();

                if (datosKmers.empty()) {
                    cout << "[ADVERTENCIA] No se encontraron k-mers para procesar.\n";
                    break;
                }

                // Construcción del Max-Heap e impresión del Top 10 de patrones más repetidos
                MaxHeapKmers maxHeap;
                maxHeap.buildHeap(datosKmers);
                maxHeap.printTopN(10);

                break;
            }
            case 5: { // Opción 5: Análisis de transiciones de 2-mers representado en un Grafo Ponderado
                cout << "\n--- GRAFO PONDERADO DE TRANSICIONES DE 2-MERS ---\n";
                if (listaSecuencias.empty()) {
                    cout << "[ERROR] No hay secuencias cargadas.\n";
                    break;
                }

                // Creación de grafo dirigido de 16 vértices (todas las combinaciones posibles de dinucleótidos)
                Grafo grafo2Mers(16, true);
                string seq = listaSecuencias[0]->sequence;

                // Estructura auxiliar local para almacenar el dinucleótido y su posición basada en 1
                struct PosKmer {
                    string kmer;
                    int posicion;
                };

                vector<PosKmer> listaPosiciones;
                int pos = 1;

                // Extracción de 2-mers avanzando paso a paso y registro de sus ubicaciones
                for (size_t i = 0; i + 1 < seq.length(); i += 2) {
                    string kmer = seq.substr(i, 2);
                    int idx = grafo2Mers.getIndex2Mer(kmer);
                    if (idx != -1) {
                        listaPosiciones.push_back({kmer, pos});
                    }
                    pos++;
                }

                double beta = 1.0; // Parámetro de atenuación para el cálculo de distancias

                // Cálculo de pesos basados en la distancia relativa (j - i) entre pares de 2-mers
                for (size_t a = 0; a < listaPosiciones.size(); ++a) {
                    for (size_t b = a + 1; b < listaPosiciones.size(); ++b) {
                        int u = grafo2Mers.getIndex2Mer(listaPosiciones[a].kmer);
                        int v = grafo2Mers.getIndex2Mer(listaPosiciones[b].kmer);
                        int i = listaPosiciones[a].posicion;
                        int j = listaPosiciones[b].posicion;

                        if (u != -1 && v != -1 && j > i) {
                            // Aplicación de la fórmula de peso basada en la inversa de la distancia
                            double peso = 1.0 / pow((j - i), beta);
                            grafo2Mers.agregarLado(u, v, peso);
                        }
                    }
                }

                cout << "Se ha construido el grafo ponderado con exito.\n";
                // Impresión de la estructura y adyacencias del grafo resultante
                grafo2Mers.imprimir();
                break;
            }
            case 6: // Opción 6: Salir de la aplicación
                cout << "\nSaliendo del sistema de analisis genomico...\n";
                break;
            default:
                cout << "\n[ADVERTENCIA] Opcion no valida. Elija un numero entre 1 y 6.\n";
        }
    } while(opcion != 6);

    // Liberación segura de la memoria dinámica reservada para las secuencias de ADN
    for (DNASequence* seq : listaSecuencias) {
        delete seq;
    }
    listaSecuencias.clear();

    return 0;
}
