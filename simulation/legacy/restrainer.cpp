#include<fstream> // std::ifstream, std::ofstream
#include<iostream> // std::cout, std::cin
#include<iomanip> // std::setprecision
#include<cmath> // sqrt
#include<vector> // std::vector<type>
#include<string> // std::string


//
std::string infile = "model_building/4tubes_10stack.pdb";
std::string outfile = "test_cpp.dat";
//
bool rings = true;
bool coms = true;
bool tubes = true;
bool torsion = true;
bool ifvari = true;
// Los tubos están colocados:
// 1 2 >>> 1 2 3 4
// 4 3
int Ntubes = 4; // Número de nanotubos
int Natom = 40; // Número de átomos en un anillo
int Nring = 10; // Número de anillos en un nanotubo
int Nca = 80; // Número de carbonos alfas de un nanotubo
bool box = false; // El pdb empieza con una caja
bool vmd = false; // El pdb ha sigo generado por vmd
//
int Ntotal = Nring*Natom; // Número total de átomos en un nanotubo
// connect1[i] está conectado con connect2[i]
const int Nconnect = 4;
int connect1[Nconnect] = {0, 1, 2, 3};//, 0, 1};
int connect2[Nconnect] = {1, 2, 3, 0};//, 2, 3};
//
int nstep2 = 10000;
//
float r3 = 19.0; // 25 o 20-19
float r3a = 20.0; // 20-19
float rk3 = 200.0; // 200 o 0 = rk2
float rk3a = rk3; // 200 o 0 = rk2a
//


int main() {

    // Leer el fichero pdb con los residuos a enlazar

    std::ifstream pdb(infile.c_str()); // Archivo .pdb
    std::string atom; // Elemento

    std::vector<int> ica; // Número del carbono alfa

    std::string sjunk; double djunk; // Para los datos que queremos ignorar
    if (box) {pdb >> sjunk >> djunk >> djunk >> djunk >> djunk >> djunk >> djunk >> sjunk >> djunk >> djunk;}
    for (int n=0; n<Nring; n++) {
        for (int i=0; i<Natom; i++) {
	    int index = Natom*n + i;

            pdb >> sjunk >> djunk >> atom >> sjunk;
	    if (vmd) {pdb >> sjunk;}
	    pdb >> djunk >> djunk >> djunk >> djunk >> djunk >> djunk;
	    if (vmd) {pdb >> sjunk;}
	    if (atom == "CA") {ica.push_back(index+1);}
	}
	if (not vmd) {pdb >> sjunk;}
    }
    
    pdb.close();


    // Escribir un el fichero de restricciones
    
    std::ofstream rst(outfile.c_str());

    if (rings) {
    rst << "#####################################\n";
    rst << "####            Rings           #####\n";
    rst << "#####################################\n";
    
    for (int n=0; n<Ntubes; n++) {
        for (int anillo=1; anillo<Nring; anillo++) { // Bucle sobre los 5 "enlaces" entre anillos del 1er tubo
	    rst << "&rst iresid = 0,\n";
	    rst << "     iat = -1, -1,\n";
	    rst << "     igr1 = ";
	    for (int i=0; i<8; i++) { // Bucle sobre los carbonos A del anterior anillo
	        rst << ica[8*(anillo-1)+i] + n*Ntotal << ",";
	    }
	    rst << "\n     igr2 = ";
	    for (int i=0; i<8; i++) { // Bucle sobre los carbonos A del anillo
	        rst << ica[8*anillo+i] + n*Ntotal << ",";
	    }
	    rst << "\n     ifvari = 0,\n";
	    rst << "     r1 = 0, r2 = 0.0, r3 = 5.0, r4 = 9999,\n";
	    rst << "     rk2 = " << rk3 << ", rk3 = " << rk3 << ",\n";
	    rst << "&end\n";
	    rst << "#####################################\n";
	}
    }
    }


    if (coms) {
    rst << "#####################################\n";
    rst << "####            COMs            #####\n";
    rst << "#####################################\n";
    
    for (int n=0; n<Nconnect; n++) {
        rst << "&rst iresid = 0,\n";
	rst << "     iat = -1, -1,\n";
        rst << "     igr1 = ";
	for (int anillo=0; anillo<Nring; anillo++) {
	    for (int i=0; i<8; i++) { // Carbonos A del 1er NT
	        rst << ica[8*anillo+i] + connect1[n]*Ntotal << ",";
	    }
	}
        rst << "\n     igr2 = ";
	for (int anillo=0; anillo<Nring; anillo++) {
	    for (int i=0; i<8; i++) { // Carbonos A del 2o NT
	        rst << ica[8*anillo+i] + connect2[n]*Ntotal << ",";
	    }
	}
	if (not ifvari) {
	rst << "\n     ifvari = 0,\n";
	rst << "     r1 = 0, r2 = 0.0, r3 = " << r3 << ", r4 = 9999,\n";
	rst << "     rk2 = " << rk3 << ", rk3 = " << rk3 << ",\n";
	rst << "&end\n";
	rst << "#####################################\n";
	}
	if (ifvari) {
	rst << "\n     ifvari = 1, nstep1 = 0, nstep2 = " << nstep2 << ",\n";
	rst << "     r1 = 0, r2 = 0.0, r3 = " << r3 << ", r4 = 9999,\n";
        rst << "     r1a = 0, r2a = 0.0, r3a = " << r3a << ", r4a = 9999,\n";
	rst << "     rk2 = " << rk3 << ", rk3 = " << rk3 << ",\n";
	rst << "     rk2a = " << rk3a << ", rk3a = " << rk3a << ",\n";
	rst << "&end\n";
	rst << "#####################################\n";
	}
    }
    }
    

    if (tubes) {
    rst << "#####################################\n";
    rst << "####         NTube-NTube        #####\n";
    rst << "#####################################\n";
    
    for (int n=0; n<Nconnect; n++) {
        for (int anillo=0; anillo<Nring; anillo+=Nring-1) {
            rst << "&rst iresid = 0,\n";
	    rst << "     iat = -1, -1,\n";
	    rst << "     igr1 = ";
	    for (int i=0; i<8; i++) { // Carbonos A del 1er NT
	        rst << ica[8*anillo+i] + connect1[n]*Ntotal << ",";
	    }
	    rst << "\n     igr2 = ";
	    for (int i=0; i<8; i++) { // Carbonos A del 2o NT
	        rst << ica[8*anillo+i] + connect2[n]*Ntotal << ",";
	    }
	    if (not ifvari) {
	    rst << "\n     ifvari = 0,\n";
	    rst << "     r1 = 0, r2 = 0.0, r3 = " << r3 << ", r4 = 9999,\n";
	    rst << "     rk2 = " << rk3 << ", rk3 = " << rk3 << ",\n";
	    rst << "&end\n";
	    rst << "#####################################\n";
	    }
	    if (ifvari) {
	    rst << "\n     ifvari = 1, nstep1 = 0, nstep2 = " << nstep2 << ",\n";
	    rst << "     r1 = 0, r2 = 0.0, r3 = " << r3 << ", r4 = 9999,\n";
	    rst << "     r1a = 0, r2a = 0.0, r3a = " << r3a << ", r4a = 9999,\n";
	    rst << "     rk2 = " << rk3 << ", rk3 = " << rk3 << ",\n";
	    rst << "     rk2a = " << rk3a << ", rk3a = " << rk3a << ",\n";
	    rst << "&end\n";
	    rst << "#####################################\n";
	    }
	}
    }
    }
    

    if (torsion) {
    rst << "#####################################\n";
    rst << "####           Torsion          #####\n";
    rst << "#####################################\n";
    
    for (int n=0; n<Nconnect; n++) {
        rst << "&rst iresid = 0,\n";
	rst << "     iat = ";
	rst << ica[8*0+0] + connect1[n]*Ntotal << ", "; //
	rst << ica[8*0+0] + connect2[n]*Ntotal << ", "; //
	rst << ica[8*(Nring-1)+0] + connect2[n]*Ntotal << ", "; //
	rst << ica[8*(Nring-1)+0] + connect1[n]*Ntotal << ","; //
	rst << "\n     ifvari = 0,\n";
	rst << "     r1 = -180, r2 = -5, r3 = 5, r4 = 180,\n";
	rst << "     rk2 = " << rk3 << ", rk3 = " << rk3 << ",\n";
	rst << "&end\n";
	rst << "#####################################\n";
    }
    }

    rst.close();

    return 0;
}
