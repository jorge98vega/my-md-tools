#include<fstream> // std::ifstream, std::ofstream
#include<iostream> // std::cin, std::cout
#include<iomanip> // std::setprecision
#include<cmath> // sqrt, cos, sin
#include<vector> // std::vector<type>
#include<string> // std::string


// Datos del pdb a clonar
int Nring = 40; // Número de anillos
int Natom = 166; // Número de átomos en cada anillo
int Nres = 8; // Número de residuos en cada anillo
bool box = true; // El pdb empieza con una caja
bool vmd = true; // El pdb ha sigo generado por vmd
//
int Ntotal = Nring*Natom; // Número total de átomos
//
std::string infile = "dry-and-reWAT/4t10s_run01_dry.pdb";
std::string outfile = "dry-and-reWAT/4t10s_run01_dry_clean.pdb";
//


void clon(std::string outfile, std::vector<std::string> atoms, std::vector<std::string> res, std::vector<int> ires, std::vector<double> pos, double com[3], double desplazo[3], double giro, int counter) {

    std::ofstream file(outfile.c_str(), std::ios_base::app);

    for (int n=0; n<Nring; n++) {
        for (int i=0; i<Natom; i++) {
	    int index = Natom*n + i;
	    int iat = counter*Ntotal + index + 1;
	    int nres = counter*Nring*Nres + ires[index];
	    
	    file << "ATOM  ";
	    if (iat < 10) {file << "    ";}
	    else if (iat < 100) {file << "   ";}
	    else if (iat < 1000) {file << "  ";}
	    else if (iat < 10000) {file << " ";}
	    file << iat << "  " << atoms[index];
	    if (atoms[index].length() < 2) {file << "  ";}
	    else if (atoms[index].length() < 3) {file << " ";}
	    file << " " << res[index] << "   ";
            if (nres < 10) {file << "  ";}
	    else if (nres < 100) {file << " ";}
	    file << nres << "     ";
	    double newpos[3] = {0.0, 0.0, 0.0};
	    newpos[0] = cos(giro)*(pos[3*index+0] - com[0]) - sin(giro)*(pos[3*index+1] - com[1]);
	    newpos[1] = sin(giro)*(pos[3*index+0] - com[0]) + cos(giro)*(pos[3*index+1] - com[1]);
	    newpos[2] = pos[3*index+2] - com[2];
	    for (int dim=0; dim<3; dim++) {
	        double newcoord = newpos[dim] + com[dim] + desplazo[dim];
		if (newcoord >= -10.0 and newcoord < 100.0) {file << " ";}
	        if (newcoord >= 0.0 and newcoord < 10.0) {file << " ";}
	        file << std::fixed << std::setprecision(3) << newcoord << " ";
            }
	    file << "0.00  0.00\n";
	}
	file << "TER   \n";
    }

    file.close();
}


int main() {
 
    // Leer el fichero del pdb
    
    std::ifstream pdb(infile.c_str()); // Archivo .pdb
    std::vector<std::string> atoms; // Elemento
    std::vector<std::string> res; // Residuo
    std::vector<int> ires; // Número de residuo
    std::vector<double> pos; // Posición del átomo
    double com[3] = {0.0, 0.0, 0.0};
    atoms.resize(Ntotal);
    res.resize(Ntotal);
    ires.resize(Ntotal);
    pos.resize(3*Ntotal);

    std::ofstream icas("CAs.dat");

    std::string sjunk; double djunk; // Para los datos que queremos ignorar
    if (box) {pdb >> sjunk >> djunk >> djunk >> djunk >> djunk >> djunk >> djunk >> sjunk >> djunk >> djunk;}
    for (int n=0; n<Nring; n++) {
        for (int i=0; i<Natom; i++) {
	    int index = Natom*n + i;
            pdb >> sjunk >> djunk >> atoms[index] >> res[index];
	    if (vmd) {pdb >> sjunk;}
	    pdb >> ires[index];
	    for (int dim=0; dim<3; dim++) {
	        pdb >> pos[3*index + dim];
		com[dim] += pos[3*index + dim];
            }
	    pdb >> djunk >> djunk;
	    if (vmd) {pdb >> sjunk;}

	    if (atoms[index] == "CA") {icas << index+1 << "\n";}
	}
	if (not vmd) {pdb >> sjunk;}
    }
    
    pdb.close();
    
    for (int dim=0; dim<3; dim++) {
        com[dim] = com[dim]/(1.0*Ntotal);
    }
    
    
    // Crear el nuevo pdb con los tubos clonados (desplazados y girados)
    
    std::ofstream newpdb(outfile.c_str());
    newpdb << "\n";
    newpdb.close();
    
    double desplazo[3] = {0.0, 0.0, 0.0}; // Desplazamiento de todos los átomos
    double giro = 0.0; // Giro en torno al eje z del tubo, en radianes

    std::ifstream input("clon.in");
    int Nclon;
    input >> Nclon;
    
    for (int n=0; n<Nclon; n++) {
      input >> desplazo [0] >> desplazo[1] >> desplazo[2] >> giro;
      clon(outfile, atoms, res, ires, pos, com, desplazo, giro, n);
    }

    input.close();

    std::ofstream newpdb2(outfile.c_str(), std::ios_base::app);
    newpdb2 << "END   \n";
    newpdb2.close();

    return 0;
}
