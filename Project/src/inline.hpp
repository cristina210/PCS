#ifndef Inline_H
#define Inline_H

#include "Eigen/Eigen"
#include "FracturesTracesPolygons.hpp"
#include <iostream>


using namespace std;
using namespace Eigen;


// Funzione per implementare il prodotto vettoriale
inline Vector3d vec_product(Vector3d& v1, Vector3d& v2){

    Vector3d v = {};
    v[0] = v1[1]*v2[2] - v2[1]*v1[2];
    v[1] = - v1[0]*v2[2] + v1[2]*v2[0];
    v[2] = v1[0]*v2[1] - v2[0]*v1[1];

    return v;
}


// Funzione per calcolare la distanza euclidea senza usare .norm()
inline double euclidean_distance(const Vector3d& a, const Vector3d& b) {
    return sqrt((a(0) - b(0)) * (a(0) - b(0)) +
                (a(1) - b(1)) * (a(1) - b(1)) +
                (a(2) - b(2)) * (a(2) - b(2)));
}

// Funzione per aggiornare le strutture per le tracce
inline void inserimento_map(int pass, unsigned int idpar, GeometryLibrary::Traces& trace) {
    if (pass == 0) {
        // Assicurati che il vettore interno esista
        if (idpar >= trace.TraceIdsPassxFracture.size()) {
            trace.TraceIdsPassxFracture.resize(idpar + 1);
        }
        // Inserisci l'ID della traccia nel vector TraceIdsPassxFracture nella posizione idpar
        trace.TraceIdsPassxFracture[idpar].push_back(trace.numTraces - 1);

    } else {
        // Assicurati che il vettore interno esista
        if (idpar >= trace.TraceIdsNoPassxFracture.size()) {
            trace.TraceIdsNoPassxFracture.resize(idpar + 1);
        }
        // Inserisci l'ID della traccia nel vector TraceIdsNoPassxFracture nella posizione idpar
        trace.TraceIdsNoPassxFracture[idpar].push_back(trace.numTraces - 1);
    }
}

// Funzione controllare se un punto è combinazione convessa di altri due
inline bool combinazione_convessa(Vector3d v1, Vector3d v2, Vector3d p){
    double alpha = 0;
    for (int i=0;i<3;i++) //presupponiamo che le fratture non siano degeneri
    {
        if (abs(v2[i]-v1[i])>tolDefault)
        {
            alpha = (p[i] - v1[i])/(v2[i] - v1[i]);
            break;
        }
    }
    if (alpha <= 1 + tolDefault && alpha >= -tolDefault) // se alpha < 0 o alpha > 1 (in realtà sbordo un pò a sx e dx dell'intervallo)
    {
        return true;
    }
    return false;


}

// Funzione usata per la memorizzazione di Cell0D in parte 2 per controllare se un punto è già esistente nella struttura dove si vuole memorizzare
///Controllo se Punto0 è presente in VettoreCoordianteIn0D --> return true se l'inserimento non è ancora avvenuto e false se è già avvenuto
inline bool checkInserimento(const Vector3d Punto0, const vector<Vector3d>& VettoreCoordinateIn0D){
    // Vettore Coordinate lo passiamo in referenza perchè potrebbe essere "grande" e quindi sarebbe oneroso fare una copia
    bool presente = false;
    for (unsigned int i = 0; i<VettoreCoordinateIn0D.size();i++)
    {
        Vector3d vettorino = VettoreCoordinateIn0D[i];//punto a 3 coordinate

        if (abs(vettorino[0] - Punto0[0])<tolDefault)//se la prima coord è uguale
        {
            if (abs(vettorino[1] - Punto0[1])<tolDefault)//se la seconda coordinata è uguale
            {
                if (abs(vettorino[2] - Punto0[2])<tolDefault) //se la terza coordinata è uguale
                {
                    presente = true;
                }
            }
        }
    }
    return presente;
}

// Funzione usata per l'aggiornamento di Cell0D nella memorizzazione di Cell0D in parte 2.
// aggiorna il numero di Cell0D,aggiunge id, coordinate e marker nel dizionario dei marker
// aggiunge una matrice matrice che corrisponde al nuovo id inserito in SequenzeXpunto
inline void addAndPrintPoint(GeometryLibrary::Polygons& sottoPoligono, map<unsigned int, list<unsigned int>>& markerDiz, const Vector3d& coordinates, unsigned int markerKey) {
    unsigned int NumPuntiFinora = sottoPoligono.NumberCell0D;
    sottoPoligono.NumberCell0D = NumPuntiFinora + 1;

    sottoPoligono.Cell0DId.push_back(NumPuntiFinora);
    sottoPoligono.Cell0DCoordinates.push_back(coordinates);

    cout << "++ Punto " << sottoPoligono.Cell0DId[NumPuntiFinora] << ": " << sottoPoligono.Cell0DCoordinates[NumPuntiFinora].transpose();
    cout << " con Marker: ";

    MatrixXd M(0, 0); // matrice vuota di dimensione
    sottoPoligono.SequenzeXpunto.push_back(M);
    markerDiz[markerKey].push_back(NumPuntiFinora); // marker con chiave markerKey
}


#endif
