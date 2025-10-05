# Progetto Discrete Fracture Network (DFN)

## Panoramica
Questo progetto implementa un'analisi di un **Discrete Fracture Network** in C++.  
Un Discrete Fracture Network è un sistema costituito da \(N\) fratture \(F_n, n=1,\dots,N\) rappresentate da poligoni planari nello spazio 3D. Le intersezioni tra fratture, chiamate **tracce** \(T_m, m=1,\dots,M\), sono modellate come segmenti.
Le tracce possono essere **passanti** (entrambi gli estremi giacciono sul bordo della frattura) o **non-passanti** (almeno un estremo si trova all'interno della frattura).  

L'obiettivo del progetto è doppio:  

1. Identificare le tracce per ogni frattura, classificarle in passanti o non-passanti e ordinarle per lunghezza.  
2. Generare i sotto-poligoni creati quando una frattura viene tagliata dalle sue tracce.

---

## Flusso principale
- Memorizzazione di vertici (0D), lati (1D) e poligoni (2D)
- Individuazioni tracce
- Taglio delle fratture secondo l’ordine:
  1. Tracce passanti
  2. Tracce non-passanti (prolungate fino a intersecare i lati del sotto-poligono)
- Preparazione dei dati per l’esportazione in **Paraview** 

---

## Struttura del codice
Il progetto è organizzato in moduli principali:

### 1. **Fractures**
Rappresenta le fratture del DFN.

### 2. **Traces**
Rappresenta le tracce generate dalle intersezioni tra fratture.

### 3. **Polygons**
Rappresenta i sotto-poligoni generati dalle fratture dopo il taglio con le tracce.
L'Oggetto `PolygonalMesh` Per ogni frattura l'oggetto : una collezione di poligoni che definiscono la superficie della frattura.

### 4. **Utility Functions**
Il progetto implementa funzioni geometriche fondamentali:
normal_vector(): calcola il vettore normale al piano definito dai primi tre vertici della frattura.
check_sphere(): verifica rapidamente se due fratture possono intersecarsi usando una condizione di sfera circoscritta.
intersezione_piani(): calcola la retta di intersezione tra due piani.
intersezione_rette(): calcola l’intersezione tra una retta e un segmento.

---

## Approccio per identificare le tracce
In **Finding_Traces**.
Le tracce vengono individuate calcolando l’intersezione tra le coppie di fratture.
Per ogni coppia:
1. Si verifica preliminarmente la possibilità di intersezione tramite il metodo delle sfere.
2. Si calcola la retta di intersezione tra i piani delle due fratture.
3. Si trovano i punti di intersezione di tale retta con i bordi di ciascuna frattura.
4. Se per entrambe vengono individuati due punti validi, si definisce la traccia come il segmento tra questi estremi.
Ogni traccia viene classificata come passante o non passante in base alla sua posizione rispetto ai bordi delle fratture.

## Approccio per la generazione dei sottopoligoni
In **Utils_part_Two**.
Per ogni frattura, i sottopoligoni vengono determinati a partire dall’intersezione tra la frattura stessa e le tracce che la attraversano (passanti e non passanti).
Il processo prevede tre fasi principali:
1. Memorizzazione dei vertici dei sottopoligoni: vengono salvati i vertici della frattura, gli estremi delle tracce passanti e i punti di intersezione tra le tracce.
2: Generazione delle sequenze: per ciascun vertice si costruisce una sequenza di bit che ne descrive la posizione relativa rispetto a ogni traccia (a destra, a sinistra o sulla traccia).
L’assegnazione del valore 0 o 1 avviene calcolando il prodotto vettoriale tra la direzione della traccia e il vettore che unisce un punto di riferimento della traccia al punto considerato:
il segno del prodotto identifica il semipiano di appartenenza.
3. Ricostruzione dei sottopoligoni: punti che condividono la stessa sequenza appartengono allo stesso sottopoligono, consentendo di suddividere la frattura in regioni distinte delimitate dalle tracce.


