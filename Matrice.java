import java.util.Scanner;
public class Matrice {
	double[][] matrice;
	int ligne;
	int colonne;
	public Matrice (int ligne, int colonne) {
        this.ligne = ligne;
        this.colonne = colonne;
        matrice = new double[ligne][colonne];
	}
	public double getValue (int ligne, int colonne){
		return matrice[ligne][colonne];
	}
	public void setValue (int ligne, int colonne, double value){
		matrice[ligne][colonne] = value;
	}
	public int getLigne() {  
        return ligne;
    }

    public int getColonne() { 
        return colonne;
    }
    public void afficherMatrice() {
        System.out.println("Matrix :  ");
        for (int i = 0; i < this.ligne; i++) {
            for (int j = 0; j < this.colonne; j++) {
				System.out.print(matrice[i][j] + "\t");
			}
		System.out.println("\n");
		}
	}
    public Matrice creerMatrice() {
    Matrice creer = new Matrice(ligne,colonne);
    Scanner scanner = new Scanner(System.in);
        for (int i = 0; i< this.ligne;i++) {
            for ( int j  = 0 ; j < this.colonne; j ++) {
                System.out.println("Element number  [" + i + ", " + j + "]: ");
                creer.setValue(i, j, scanner.nextInt());
            }
        }
        return creer;
    }
    
    public boolean testCarre() {
		if(getLigne() == getColonne()){    
			return true;
		}
		else {
			return false;
		}
	}

    public Matrice transposeMatrice() {
		Matrice transposeeMatrice = new Matrice(this.colonne,this.ligne);                                            
        for (int i =0; i<this.ligne; i++) {
            for (int j =0; j<this.colonne; j ++) {
                transposeeMatrice.setValue(j, i, this.getValue(i,j));                    
            }
        }
        return transposeeMatrice;
    }

    public Matrice addition(Matrice m2){                    
		int ligne1 = getLigne();
		int colonne1 = getColonne();
		Matrice result = new Matrice(ligne1, colonne1);
		for (int i = 0; i < ligne1; i++) {
			for (int j = 0; j < colonne1; j++) {
				result.setValue(i, j, (int) (this.getValue(i, j) + m2.getValue(i, j)));
			}
		}
		return result;
	}
	public Matrice multiplication(Matrice m2) {
		int ligne2 = this.getLigne();                                   
		int colonne2 = m2.getColonne();
		Matrice result = new Matrice(ligne2, colonne2);
		for (int i = 0; i < ligne2; i++) {
			for (int j = 0; j < colonne2; j++) {
				for (int k = 0; k < m2.getLigne(); k++) {
					result.setValue(i, j, (int) (result.getValue(i, j) + this.getValue(i, k) * m2.getValue(k, j)));
				}
			}
		}
		return result;
	}
	
	// transferer matrice --> array
	public double[][] transfererAArray(Matrice m1) {
		double a1 [][] = new double [this.ligne][this.colonne];
		for (int i = 0; i <this.ligne; i++ ) {
			for (int j=0; j<this.colonne;j++){
				a1[i][j] =  m1.getValue(i,j);
			}
		}
    return  a1;
	}
	public double detCofacteur(double matrice[][]){
		double temporary[][] ;
		double result = 0;   
		if (matrice.length == 1) {
			result = matrice[0][0];
			return result;
		}
		
		if (matrice.length == 2) {
			result = ((matrice[0][0] * matrice[1][1]) - (matrice[0][1] * matrice[1][0]));
			return (double) (result);
		}
		
		for (int i = 0; i < matrice[0].length; i++) {
			temporary = new double[matrice.length - 1][matrice[0].length - 1];
			for (int j = 1; j < matrice.length; j++) {
				for (int k = 0; k < matrice[0].length; k++) {
					if (k < i) {
						temporary[j - 1][k] = matrice[j][k];
					} else if (k > i) {
						temporary[j - 1][k - 1] = matrice[j][k];
					}
				}
			}
				
			result += matrice[0][i] * Math.pow (-1, (double) i)*detCofacteur(temporary);
		}
		return (double) (result);
    }
    
    // for methode detGauss
	public Matrice ligneEchangerGauss(int ligneX, int ligneY) {   
		double[] echange = new double[getLigne()];
		for (int i = 0; i<this.ligne; i++) {         
			echange[i] = this.getValue(ligneX, i);
			this.setValue(ligneX, i, -this.getValue(ligneY, i));
			this.setValue(ligneY, i, echange[i]);
		}
		return this;
	}
	public double detGauss() {
		if (this.testCarre() == true) {
			int n = getLigne();
			// creer une pseudo matrice
			Matrice pseudo = new Matrice(n, n);
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < n; j++) {
					pseudo.setValue(i,j,this.getValue(i,j));
				}
			}	 
			int indicePivot = 0;
			double[] coefEliminer = new double[n];
			double det = 1;
			// elimination  Gauss, create 1 matrix trianglulaire
			while (indicePivot < n) {
				// trouver indice pivot != 0
				if (pseudo.getValue(indicePivot,indicePivot) == 0) {
					for (int i = indicePivot; i < n; i++) {   
						if (pseudo.getValue(indicePivot, i) != 0) {
							pseudo.ligneEchangerGauss(indicePivot, i);
							break;
						}
						if (i == n-1) {
							return det = 0;
						}
					}
					// si indice pivot != 0 n'existe pas, determinant = 0	
				}
				// finding coefficient of elimination
				for (int i = indicePivot + 1; i < n; i++) {            
					coefEliminer[i] = pseudo.getValue(i, indicePivot)/pseudo.getValue(indicePivot,indicePivot);
				}
				// elimination
				for (int i = indicePivot + 1; i < n; i++) {            
					for (int j = indicePivot; j < n; j++) {
						pseudo.setValue(i, j, pseudo.getValue(i,j) - coefEliminer[i]*pseudo.getValue(indicePivot, j) );
					}
				}
				indicePivot++;
			}
			// --> now, matrix pseudo is a matrice triangulaire 
			// calculate determinant
			for (int i = 0; i < n; i++) {
				det *= pseudo.getValue(i,i);
			}
			return det;
		}
		else {
			System.out.println("Matrix isn't carree doesn't have determinant");
			return 0;
		}
	}
	
	// pour method inverse()
	public Matrice ligneEchangerInverse(int ligneX, int ligneY) {   
		double[] echange = new double[getLigne()];
		for (int i = 0; i<this.ligne; i++) {         
			echange[i] = this.getValue(ligneX, i);
			this.setValue(ligneX, i, this.getValue(ligneY, i));
			this.setValue(ligneY, i, echange[i]);
		}
		return this;
	}
	public Matrice identiteMatrice(int degree) {
		Matrice iD = new Matrice(degree,degree);
		for (int i = 0; i<degree; i++) {
			iD.setValue(i,i,1);
		}
		return iD;
	}
    public Matrice arrondir(int precision) {
		int scale = (int) Math.pow(10,precision);
		for (int i = 0; i < this.ligne; i++) {
			for (int j = 0; j < this.colonne; j++) {
				// -0.0 --> 0.0
				if (this.getValue(i,j) == -0) 
					this.setValue(i,j,0);
				this.setValue(i, j, (double) Math.round(this.getValue(i,j)*scale)/scale);
			}
		}
		return this;
	}
	public Matrice inverse() {
		int n = this.ligne; 
		Matrice pseudo = new Matrice(n, n);
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				pseudo.setValue(i,j,this.getValue(i,j));
			}
		}
		Matrice iD = identiteMatrice(n); 
		int indicePivot = 0;
		double coefDiviser = 0;
		double[] coefEliminer = new double[n];
		while (indicePivot < n) {
			if (pseudo.getValue(indicePivot,indicePivot) == 0) {
				for (int i = indicePivot; i < n; i++) {   
					if (pseudo.getValue(i,indicePivot) != 0) {
						pseudo.ligneEchangerInverse(indicePivot, i);
						iD.ligneEchangerInverse(indicePivot, i);
						break;
					}
				}
			}
			// divise --> matrice[i][i] = 1
			if (pseudo.getValue(indicePivot,indicePivot) != 1) {   
				coefDiviser = pseudo.getValue(indicePivot, indicePivot);
				for (int i = 0; i < n; i++) {
					pseudo.setValue( indicePivot, i, (pseudo.getValue(indicePivot,i)/coefDiviser ));
					iD.setValue( indicePivot, i, (iD.getValue(indicePivot,i)/coefDiviser ));
				}
			}
			// finding coeficient
			for (int i = 0; i < n; i++) {           
				if (i != indicePivot) { 
					coefEliminer[i] = pseudo.getValue(i, indicePivot);
				}
			}
			// elimination
			for (int i = 0; i < n; i++) {           
				if (i != indicePivot) {
					for (int j = 0; j < n; j++) {
						pseudo.setValue(i, j, pseudo.getValue(i,j) - coefEliminer[i]*pseudo.getValue(indicePivot, j) );
						iD.setValue(i, j, iD.getValue(i,j) - coefEliminer[i]*iD.getValue(indicePivot, j) );
					}
				}
			}
			indicePivot++;
		}
		iD.arrondir(2);
		return iD;
	}
	
}
