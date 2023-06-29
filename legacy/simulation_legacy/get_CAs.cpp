// (Este programa ya no es necesario)
// Uso:
// g++ get_CAs.cpp -o get_CAs
// ./get_CAs


#include<fstream> // std::ifstream, std::ofstream
#include<iostream> // std::cin, std::cout
#include<iomanip> // std::setprecision
#include<cmath> // sqrt, cos, sin
#include<vector> // std::vector<type>


// Datos del pdb
int Nring = 6; // Número de anillos en el pdb
int Natom = 164; // Número de átomos en cada anillo
bool box = true; // El pdb empieza con una caja
bool vmd = false; // El pdb ha sigo generado por vmd
//


int main() {
 
    // Leer el fichero del pdb
    std::ifstream pdb("4tubes_run03.pdb"); // Archivo .pdb

    // Escribir el archivo con los índices de los carbonos alfa
    std::ofstream icas("CAs.dat");
    std::string atom; // Elemento del átomo
    std::string sjunk; double djunk; // Para los datos que queremos ignorar

    if (box) {pdb >> sjunk >> djunk >> djunk >> djunk >> djunk >> djunk >> djunk >> sjunk >> djunk >> djunk;}
    for (int n=0; n<Nring; n++) {
        for (int i=0; i<Natom; i++) {
	    int index = Natom*n + i;

            pdb >> sjunk >> djunk >> atom >> sjunk;
	    if (vmd) {pdb >> sjunk;}
	    pdb >> djunk >> djunk >> djunk >> djunk >> djunk >> djunk;
	    if (vmd) {pdb >> sjunk;}

	    if (atom == "CA") {icas << index+1 << "\n";}
	}
	if (not vmd) {pdb >> sjunk;}
    }
    
    pdb.close();
    
    return 0;
}
