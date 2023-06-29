#include<fstream> // std::ifstream, std::ofstream
#include<iomanip> // std::setprecision
#include<cmath> // sqrt
#include<vector> // std::vector<type>


int main() {

    // Leer el fichero con los residuos a enlazar
  
    std::ifstream hbonds("Hbonds_1.txt"); // Archivo a leer
    int Nres = 40; // Número de residuos en el fichero
    std::vector<int> ires; // Número del residuo
    ires.resize(Nres);

    for (int i=0; i<Nres; i++) {
        hbonds >> ires[i];
    }

    hbonds.close();


    // Escribir un el fichero de restricciones
    
    std::ofstream rst("Hbonds_rst.dat");

    rst << "#####################################\n";
    rst << "####             HB             #####\n";
    rst << "#####################################\n";

    for (int n=0; n<20; n++) { // Bucle sobre los 20 enlaces
	rst << "&rst iresid = 1,\n";
        rst << "     iat = " << ires[2*n] << ", " << ires[2*n+1] << ",\n";
	rst << "     atnam(1) = 'O', atnam(2) = 'N',\n";
	rst << "     ifvari = 0,\n";
	//rst << "     ifvari = 1, nstep1 = 0, nstep2 = 100000,\n";
	rst << "     r1 = 0, r2 = 0.0, r3 = 3.0, r4 = 9999\n";
	//rst << "     r1a = 0, r2a = 0.0, r3a = 3.0, r4a = 9999,\n";
	rst << "     rk2 = 200, rk3 = 200,\n";
	//rst << "     rk2a = 0, rk3a = 0,\n";
	rst << "&end\n";
	rst << "&rst iresid = 1,\n";
        rst << "     iat = " << ires[2*n] << ", " << ires[2*n+1] << ",\n";
	rst << "     atnam(1) = 'N', atnam(2) = 'O',\n";
	rst << "     ifvari = 0,\n";
	//rst << "     ifvari = 1, nstep1 = 0, nstep2 = 100000,\n";
	rst << "     r1 = 0, r2 = 0.0, r3 = 3.0, r4 = 9999\n";
	//rst << "     r1a = 0, r2a = 0.0, r3a = 3.0, r4a = 9999,\n";
	rst << "     rk2 = 200, rk3 = 200,\n";
	//rst << "     rk2a = 0, rk3a = 0,\n";
	rst << "&end\n";
	rst << "#####################################\n";
    }

    for (int n=0; n<20; n++) { // Bucle sobre los 20 enlaces de H del 2o tubo
        rst << "&rst iresid = 1,\n";
        rst << "     iat = " << ires[2*n]+48 << ", " << ires[2*n+1]+48 << ",\n";
	rst << "     atnam(1) = 'O', atnam(2) = 'N',\n";
	//rst << "     ifvari = 0,\n";
	rst << "     ifvari = 1, nstep1 = 0, nstep2 = 100000,\n";
	rst << "     r1 = 0, r2 = 0.0, r3 = 3.0, r4 = 9999,\n";
	rst << "     r1a = 0, r2a = 0.0, r3a = 3.0, r4a = 9999,\n";
	rst << "     rk2 = 200, rk3 = 200,\n"; // 200 si ifvari = 1, 0 si no
	rst << "     rk2a = 0, rk3a = 0,\n";
	rst << "&end\n";
	rst << "&rst iresid = 1,\n";
	rst << "     iat = " << ires[2*n]+48 << ", " << ires[2*n+1]+48 << ",\n";
	rst << "     atnam(1) = 'N', atnam(2) = 'O',\n";
	//rst << "     ifvari = 0,\n";
	rst << "     ifvari = 1, nstep1 = 0, nstep2 = 100000,\n";
	rst << "     r1 = 0, r2 = 0.0, r3 = 3.0, r4 = 9999,\n";
	rst << "     r1a = 0, r2a = 0.0, r3a = 3.0, r4a = 9999,\n";
	rst << "     rk2 = 200, rk3 = 200,\n"; // 200 si ifvari = 1, 0 si no
	rst << "     rk2a = 0, rk3a = 0,\n";
	rst << "&end\n";
	rst << "#####################################\n";
    }

    rst.close();

    return 0;
}
