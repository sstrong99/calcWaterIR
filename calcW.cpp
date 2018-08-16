#include "calcW.h"

CalcW::CalcW(const int model,const int natoms,float avef) : avef(avef) {
  float qO;
  if (model==0) { //SPCE
    aPerM=3;
    qO=-0.8476f;
    mInc=0;
    mymap = new MapAS2008();
    h2x = new SPCEind();
  } else if (model==1 || model==2) { //TIP4P class
    aPerM=4;
    qO=-1.04f;
    mInc=aPerM-1;
    mymap = new MapGruenbaum2013();
    h2x = new TIP4Pind();
  } else {
    printf("ERROR: invalid model choice\n");
    exit(EXIT_FAILURE);
  }

  charges = new float[aPerM];
  if (model==0) {
    charges[0]=qO;
    charges[1]=-charges[0]/2;
  } else {
    charges[0]=0.0f;
    charges[3]=qO;
    charges[1]=-charges[3]/2;
  }
  charges[2]=charges[1];

  moveMfrac=0.0;
  if (model==2) //for TIP4P/2005 and E3B3
    moveMfrac=0.15/0.1546;

  if (natoms % aPerM) {
    printf("ERROR: the number of atoms is not divisible by %d\n",aPerM);
    exit(EXIT_FAILURE);
  }
  nO=natoms/aPerM;
#ifdef DEBUG
  nO=30;
#endif
  nH=nO*2;

  E      = new float[nH];
  OH     = new rvec[nH];
  dipdip = new float[nH*nH];
  wMat   = new float[nH*nH];

  prnt.setNH(nH);
}

CalcW::~CalcW() {
  delete[] E;
  delete[] OH;
  delete[] dipdip;
  delete[] wMat;
  delete[] charges;
  delete mymap;
  delete h2x;
}

void CalcW::compute(Traj &traj, rvec *m) {
  if (moveMfrac != 0.0)
    traj.moveM(moveMfrac,aPerM);

  calcE(traj);

  //map E-field to hamiltonian elements
  mapE2W(m);
}

//This is slower than calcE, and only calculates the E-field at the H atom
//not the dipole-dipole term
void CalcW::calcE(const Traj &traj) {
  float cut2=mymap->getcut2()*A0INV*A0INV;
  float Ddip=mymap->getDdip()*A0INV;
  const rvec *x=traj.getCoords();
  rvec box;
  traj.getBox(box);
  int ii;
  float tmpcut=2*mymap->getcut()*A0INV;
  //check that box is larger than 2*cutoff
  for (ii=0; ii<DIM; ii++)
    if (box[ii]<tmpcut) {
      printf("ERROR: Box is smaller than twice the cutoff\n");
      exit(EXIT_FAILURE);
    }

  //compute OH vectors and dip locations
  rvec OHtmp;
  float d2,d;
  int hInd,jj;
  rvec *dip=new rvec[nH];  //positions of point dipoles
  for (ii=0; ii<nO; ii++) {
    for (jj=0; jj<2; jj++) { //loop through 2 H on each oxygen
      hInd=ii*aPerM + jj + 1;
      addRvec(x[hInd],x[ii*aPerM],OHtmp,-1);
      pbc(OHtmp,box);
      d2=norm2vec(OHtmp);
      d=sqrt(d2);
      multRvec(OHtmp,1.0/d);
      setRvec(OH[ii*2+jj],OHtmp);

      multRvec(OHtmp,Ddip); //OH is unit, so now has length Ddip
      addRvec(x[ii*aPerM],OHtmp,+1);
      setRvec(dip[ii*2+jj],OHtmp);
    }
  }

  //init dipdip to 0
  //TODO: this should be unecessary
  for (ii=0; ii<nH*nH; ii++)
    dipdip[ii]=0.0;

  int hi,hj,kk;
  rvec hiv,vec,dipI,OHi,tmpEi;
  float dipdiptmp;
  for (ii=0; ii<nH; ii++) { //loop over all Hs
    hi=h2x->convert(ii);

    setRvec(tmpEi,0.0);
    getRvec(x[hi],hiv);
    getRvec(OH[ii],OHi);
    getRvec(dip[ii],dipI);

    //loop through other molecules
    for (jj=0; jj<nO; jj++) {
      if ( floor(ii/2) == jj ) //skip same molecule
	continue;

      //get OH distance
      addRvec(hiv,x[jj*aPerM],vec,-1); //points from O to H
      pbc(vec,box);
      d2=norm2vec(vec);
      if (d2 < cut2) { //NOTE: cutoff is w.r.t. O position of each molecule
	d=sqrt(d2);
	multRvec(vec, charges[0]/(d*d*d) );
	addRvec(vec,tmpEi,+1);
	for (kk=1; kk<aPerM; kk++) { //loop through other atoms
	  addRvec(hiv,x[jj*aPerM+kk],vec,-1); //points from other to H
	  pbc(vec,box);
	  d2=norm2vec(vec);
	  d=sqrt(d2);

	  multRvec(vec, charges[kk]/(d*d*d));
	  addRvec(vec,tmpEi,+1);
	}
      }

      //compute dipdip coupling
      if (jj*2>ii) { //only compute upper triangular part
	for (kk=0; kk<2; kk++) { //loop over 2 Hs per molecule
	  hj=2*jj+kk;
	  addRvec(dipI,dip[hj],vec,-1);
	  pbc(vec,box);
	  d2=norm2vec(vec);
	  d=sqrt(d2);
	  multRvec(vec,1.0/d);
	  dipdiptmp=dot(OHi,OH[hj]) - 3*dot(OHi,vec)*dot(OH[hj],vec);
	  dipdiptmp/=d*d*d;
	  setNN(dipdip,dipdiptmp,ii,hj);
	}
      }
    }
    E[ii]=dot(tmpEi,OH[ii]);
  }

  delete[] dip;
}

void CalcW::calcE_wrong(const Traj &traj) {
  float cut2=mymap->getcut2()*A0INV*A0INV;
  float Ddip=mymap->getDdip()*A0INV;
  const rvec *x=traj.getCoords();
  rvec box;
  traj.getBox(box);

  rvec OHtmp,OHdip;
  float d2,d;
  int hInd,ii,jj;
  rvec *dip=new rvec[nH];  //positions of point dipoles
  rvec *Evec=new rvec[nH];
  //zero Efield and dipdip
  for (ii=0; ii<nH; ii++) {
    setRvec(Evec[ii],0.0);
    for (jj=0; jj<nH; jj++)
      setNN(dipdip,0.0,ii,jj);
  }

  //compute OH vectors and dip locations
  for (ii=0; ii<nO; ii++) {
    for (jj=0; jj<2; jj++) { //loop through 2 H on each oxygen
      hInd=ii*aPerM + jj + 1;
      addRvec(x[hInd],x[ii*aPerM],OHtmp,-1);
      pbc(OHtmp,box);
      d2=norm2vec(OHtmp);
      d=sqrt(d2);
      multRvec(OHtmp,1.0/d);
      setRvec(OH[ii*2+jj],OHtmp);

      multRvec(OHtmp,Ddip); //OH is unit, so now has length Ddip
      addRvec(x[ii*aPerM],OHtmp,OHdip,+1);
      setRvec(dip[ii*2+jj],OHdip);
    }
  }

  //compute field from HH interactions
  int hi,hj,oi,oj;
  rvec hiv,hjv,vec;
  rvec dipI,dipJ,dvec,OHi,OHj;
  float dipdiptmp;
  rvec ojv;
  for (ii=0; ii<nH-1; ii++) {
    hi=h2x->convert(ii);  //index in atom array
    getRvec(x[hi],hiv);
    getRvec(dip[ii],dipI);
    getRvec(OH[ii],OHi);
    rvec tmpEi={0.0,0.0,0.0};

    //loop through other Hs, skipping same molecule
    for (jj=ii+1+(ii+1)%2; jj<nH; jj++) {
      hj=h2x->convert(jj);  //index in atom array
      getRvec(x[hj],hjv);

      addRvec(hiv,hjv,vec,-1);
      pbc(vec,box);
      d2=norm2vec(vec);
      if (d2<cut2) {
	d=sqrt(d2);
	multRvec( vec , charges[1]/(d*d*d) );
	addRvec(vec,tmpEi,+1);
	addRvec(vec,Evec[jj],-1);
      }

      getRvec(dip[jj],dipJ);
      getRvec(OH[jj],OHj);

      addRvec(dipI,dipJ,dvec,-1);
      pbc(dvec,box);
      d2=norm2vec(dvec);
      d=sqrt(d2);
      multRvec(dvec,1.0/d);
      dipdiptmp=dot(OHi,OHj) - 3*dot(OHi,dvec)*dot(OHj,dvec);
      dipdiptmp/=d*d*d;
      setNN(dipdip,dipdiptmp,ii,jj);
      //dipdip is only accessed in the same format of loop, so don't need
      //to set both elements symetrically
      //setNN(dipdip,dipdiptmp,jj,ii);
    }

    //compute Efield from O atoms
    oi=aPerM*floor(hi/aPerM) + mInc; //O ind of this molec, to be skipped
    for (jj=0; jj<nO; jj++)
    {
      oj=jj*aPerM + mInc; //index of O in atom array
      if (oj==oi) continue;  //skip oi
      getRvec(x[oj],ojv);
      addRvec(hiv,ojv,vec,-1);
      pbc(vec,box);
      d2=norm2vec(vec);
      if (d2<cut2)
      {
	d=sqrt(d2);
	multRvec( vec , charges[mInc]/(d*d*d) );
	addRvec(vec,tmpEi,+1);
      }
    }
    addRvec(tmpEi,Evec[ii],+1);
  }

  //add last OH interaction
  ii=nH-1;
  hi=h2x->convert(ii);  //index in atom array
  getRvec(x[hi],hiv);
  rvec tmpEi={0.0,0.0,0.0};
  for (jj=0; jj<nO-1; jj++)  //can skip this O explicitely
  {
    oj=jj*aPerM+mInc; //index of O in atom array
    getRvec(x[oj],ojv);
    addRvec(hiv,ojv,vec,-1);
    pbc(vec,box);
    d2=norm2vec(vec);
    if (d2<cut2)
    {
      d=sqrt(d2);
      multRvec( vec , charges[mInc]/(d*d*d) );
      addRvec(vec,tmpEi,+1);
    }
  }
  addRvec(tmpEi,Evec[ii],+1);

  delete[] dip;

  //compute abs. mag of dot product of OH with E
  for (ii=0; ii<nH; ii++)
    E[ii]=dot(Evec[ii],OH[ii]);
  delete[] Evec;
}

float CalcW::calcAveF() {
  if (avef != 0.0)
    printf("WARNING: avef has already been set\n");

  float sum=0.0;
  for (int ii=0; ii<nH; ii++)
    sum += mymap->w(E[ii]);
  avef = sum/(float) nH;
  return avef;
}

void CalcW::mapE2W(rvec *m) {
  int ii,jj,kk;
  float wIntra,tmpw,tmpk,tmpx,tmpmu;
  float *xh  = new float[nH];
  float *mud = new float[nH];
  float *p   = new float[nH];

  //prnt.vec(E);

  //compute various map quantities and transition dipoles
  for (ii=0; ii<nH; ii++)
    {
      tmpw=mymap->w(E[ii]);
      tmpmu=mymap->mud(E[ii]);
      tmpx=mymap->x(tmpw);

      setRvec(m[ii],OH[ii]);
      multRvec(m[ii],tmpmu*tmpx);

      p[ii]=mymap->p(tmpw);
      mud[ii]=tmpmu;
      xh[ii]=tmpx;
      setNN(wMat,(tmpw-avef)*CM2PS,ii,ii);
    }

  //loop through Hs on same molecule to compute INTRAmolecluar couplings
  for (kk=0; kk<nH/2; kk++)
    {
      //get H inds on the same
      ii=kk*2;
      jj=ii+1;
      wIntra=mymap->kintra(E[ii],E[jj],xh[ii],xh[jj],p[ii],p[jj])*CM2PS;
      setNN(wMat,wIntra,ii,jj);
      setNN(wMat,wIntra,jj,ii);
    }

  //intermolecular couplings
  //loop through all H pairs, skipping those on same molecule
  for (ii=0; ii<nH-1; ii++)
    for (jj=ii+1+(ii+1)%2; jj<nH; jj++)
      {
	tmpk=xh[ii]*xh[jj]*mud[ii]*mud[jj]*getDipdip(ii,jj)*HART2CM*CM2PS;
	setNN(wMat,tmpk,ii,jj);
	setNN(wMat,tmpk,jj,ii);
      }

  delete[] xh;
  delete[] mud;
  delete[] p;

  //printf("inter,intra= %.1f\t%.1f\n",maxInterFreq(),avgIntraFreq());
  //prnt.mat(wMat);
}

float CalcW::maxInterFreq() const {
  int jj;
  int skipi;
  float max,tmp;
  float sum=0;
  float max1=0.0;
  for (int ii=0; ii<nH; ii++)
  {
    skipi=2*(ii/2)+1-(ii%2);
    max=0;
    for (jj=0; jj<nH; jj++)
    {
      if ( jj==ii || jj == skipi )
	continue;
      else
      {
	tmp=wMat[jj+ii*nH];
	if (fabs(tmp) > max)
	{
	  max=fabs(tmp);
	  max1=tmp;
	}
      }
    }
    sum+=max1;
  }
  return sum/nH/CM2PS;
}

float CalcW::avgIntraFreq() const {
  int jj;
  float sum=0;
  for (int ii=0; ii<nH; ii+=2)
  {
    jj=ii+1;
    sum+=wMat[jj+ii*nH];
  }
  return sum/(nH/2)/CM2PS;
}

void CalcW::getW(float* out) const {
  for (int ii=0; ii<nH*nH; ii++)
    out[ii]=wMat[ii];
}
