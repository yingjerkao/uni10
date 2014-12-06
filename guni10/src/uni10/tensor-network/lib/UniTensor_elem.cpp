/****************************************************************************
*  @file CMakeLists.txt
*  @license
*    Universal Tensor Network Library
*    Copyright (c) 2013-2014
*    Yun-Da Hsieh, Pochung Chen and Ying-Jer Kao
*
*    This file is part of Uni10, the Universal Tensor Network Library.
*
*    Uni10 is free software: you can redistribute it and/or modify
*    it under the terms of the GNU Lesser General Public License as published by
*    the Free Software Foundation, either version 3 of the License, or
*    (at your option) any later version.
*
*    Uni10 is distributed in the hope that it will be useful,
*    but WITHOUT ANY WARRANTY; without even the implied warranty of
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*    GNU Lesser General Public License for more details.
*
*    You should have received a copy of the GNU Lesser General Public License
*    along with Uni10.  If not, see <http://www.gnu.org/licenses/>.
*  @endlicense
*  @brief Main specification file for CMake
*  @author Ying-Jer Kao
*  @date 2014-05-06
*  @since 0.1.0
*
*****************************************************************************/
#include <uni10/tensor-network/UniTensor.h>
#include <uni10/tools/uni10_tools.h>
#include <uni10/numeric/uni10_lapack.h>
//using namespace uni10::datatype;
namespace uni10{
void UniTensor::grouping(){
	blocks.clear();
	int row_bondNum = 0;
	int col_bondNum = 0;
	RQdim = 1;
	CQdim = 1;
        bool IN_BONDS_BEFORE_OUT_BONDS = true;
	for(int i = 0; i < bonds.size(); i++){
	  if(bonds[i].type() == BD_IN){
            if(!(IN_BONDS_BEFORE_OUT_BONDS == true)){
              std::ostringstream err;
              err<<"Error in the input bond array: BD_OUT bonds must be placed after all BD_IN bonds.";
              throw std::runtime_error(exception_msg(err.str()));
            }
	    RQdim *= bonds[i].Qnums.size();
	    row_bondNum++;
	  }
	  else{
            CQdim *= bonds[i].Qnums.size();
            col_bondNum++;
            IN_BONDS_BEFORE_OUT_BONDS = false;
	  }
	}
	RBondNum = row_bondNum;
	std::map<Qnum,size_t> row_QnumMdim;
	std::vector<int> row_offs(row_bondNum, 0);
	std::map<Qnum,std::vector<int> > row_Qnum2Qidx;
	Qnum qnum;
        size_t dim;
	int boff = 0;
	std::vector<size_t>tmpRQidx2Dim(RQdim, 1);
	std::vector<size_t>tmpCQidx2Dim(CQdim, 1);
	std::vector<size_t>tmpRQidx2Off(RQdim, 0);
	std::vector<size_t>tmpCQidx2Off(CQdim, 0);
	if(row_bondNum){
		while(1){
			qnum.assign();
			dim = 1;
			for(int b = 0; b < row_bondNum; b++){
				qnum = qnum * bonds[b].Qnums[row_offs[b]];
				dim *= bonds[b].Qdegs[row_offs[b]];
			}
			if(row_QnumMdim.find(qnum) != row_QnumMdim.end()){
				tmpRQidx2Off[boff] = row_QnumMdim[qnum];
				tmpRQidx2Dim[boff] = dim;
				row_QnumMdim[qnum] += dim;
			}
			else{
				tmpRQidx2Off[boff] = 0;
				tmpRQidx2Dim[boff] = dim;
				row_QnumMdim[qnum] = dim;
			}
			row_Qnum2Qidx[qnum].push_back(boff);
			boff++;
			int bidx;
			for(bidx = row_bondNum - 1; bidx >= 0; bidx--){
				row_offs[bidx]++;
				if(row_offs[bidx] < bonds[bidx].Qnums.size())
					break;
				else
					row_offs[bidx] = 0;
			}
			if(bidx < 0)	//run over all row_bond offsets
				break;
		}
	}
	else{
		qnum.assign();
		row_QnumMdim[qnum] = 1;
		row_Qnum2Qidx[qnum].push_back(0);
	}
	std::map<Qnum,size_t> col_QnumMdim;
	std::vector<int> col_offs(col_bondNum, 0);
	std::map<Qnum,std::vector<int> > col_Qnum2Qidx;
	boff = 0;
	if(col_bondNum){
		while(1){
			qnum.assign();
			dim = 1;
			for(int b = 0; b < col_bondNum; b++){
				qnum = qnum * bonds[b + row_bondNum].Qnums[col_offs[b]];
				dim *= bonds[b + row_bondNum].Qdegs[col_offs[b]];
			}
			if(row_QnumMdim.find(qnum) != row_QnumMdim.end()){
				if(col_QnumMdim.find(qnum) != col_QnumMdim.end()){
					tmpCQidx2Off[boff] = col_QnumMdim[qnum];
					tmpCQidx2Dim[boff] = dim;
					col_QnumMdim[qnum] += dim;
				}
				else{
					tmpCQidx2Off[boff] = 0;
					tmpCQidx2Dim[boff] = dim;
					col_QnumMdim[qnum] = dim;
				}
				col_Qnum2Qidx[qnum].push_back(boff);
			}
			boff++;
			int bidx;
			for(bidx = col_bondNum - 1; bidx >= 0; bidx--){
				col_offs[bidx]++;
				if(col_offs[bidx] < bonds[bidx + row_bondNum].Qnums.size())
					break;
				else
					col_offs[bidx] = 0;
			}
			if(bidx < 0)	//run over all row_bond offsets
				break;
		}
	}
	else{
		qnum.assign();
		if(row_QnumMdim.find(qnum) != row_QnumMdim.end()){
			col_QnumMdim[qnum] = 1;
			col_Qnum2Qidx[qnum].push_back(0);
		}
	}

	std::map<Qnum,size_t>::iterator it;
	std::map<Qnum,size_t>::iterator it2;
	std::set<int> Qidx;
	Block blk;
	int qidx;
	int off = 0;
	for ( it2 = col_QnumMdim.begin() ; it2 != col_QnumMdim.end(); it2++ ){
		it = row_QnumMdim.find(it2->first);
		blk.Rnum = it->second;
		blk.Cnum = it2->second;
		blk.qnum = it2->first;
		blk.offset = off;
		off += blk.Rnum * blk.Cnum;
		blocks[it->first] = blk;
		Block* blkptr= &(blocks[it->first]);
		std::vector<int>& tmpRQidx = row_Qnum2Qidx[it->first];
		std::vector<int>& tmpCQidx = col_Qnum2Qidx[it->first];
		for(int i = 0; i < tmpRQidx.size(); i++){
			RQidx2Blk[tmpRQidx[i]] = blkptr;
			for(int j = 0; j < tmpCQidx.size(); j++){
				RQidx2Dim[tmpRQidx[i]] = tmpRQidx2Dim[tmpRQidx[i]];
				RQidx2Off[tmpRQidx[i]] = tmpRQidx2Off[tmpRQidx[i]];
				CQidx2Dim[tmpCQidx[j]] = tmpCQidx2Dim[tmpCQidx[j]];
				CQidx2Off[tmpCQidx[j]] = tmpCQidx2Off[tmpCQidx[j]];
				qidx = tmpRQidx[i] * CQdim + tmpCQidx[j];
				Qidx.insert(qidx);
			}
		}
	}
	size_t elemEnc = 0;
	for(std::map<int, size_t>::iterator itr = RQidx2Dim.begin(); itr != RQidx2Dim.end(); itr++)
		for(std::map<int, size_t>::iterator itc = CQidx2Dim.begin(); itc != CQidx2Dim.end(); itc++){
			qidx = itr->first * CQdim + itc->first;
			if(Qidx.find(qidx) != Qidx.end()){
				QidxEnc[qidx] = elemEnc;
				elemEnc += RQidx2Dim[itr->first] * CQidx2Dim[itc->first];
			}
		}
}

void UniTensor::addGate(const std::vector<_Swap>& swaps){
  try{
  if((status & HAVEBOND) == 0){
    std::ostringstream err;
    err<<"Adding swap gates to a tensor without bonds(scalar).";
    throw std::runtime_error(exception_msg(err.str()));
  }
	if((status & HAVEELEM) == 0){
    std::ostringstream err;
    err<<"Cannot add swap gates to a tensor before setting its elements.";
    throw std::runtime_error(exception_msg(err.str()));
  }
	int sign = 1;
	int bondNum = bonds.size();
	std::vector<int> Q_idxs(bondNum, 0);
	std::vector<int> Q_Bdims(bondNum, 0);
	for(int b = 0; b < bondNum; b++)
		Q_Bdims[b] = bonds[b].Qnums.size();
	int Q_off;
	int tmp;
	int RQoff, CQoff;
	size_t sB_r, sB_c;	//sub-block of a Qidx
	size_t sB_rDim, sB_cDim;	//sub-block of a Qidx
	size_t B_cDim;
	size_t E_off;
	for(std::map<int, size_t>::iterator it = QidxEnc.begin(); it != QidxEnc.end(); it++){
		Q_off = it->first;
		tmp = Q_off;
		for(int b = bondNum - 1; b >= 0; b--){
			Q_idxs[b] = tmp % Q_Bdims[b];
			tmp /= Q_Bdims[b];
		}
		RQoff = Q_off / CQdim;
		CQoff = Q_off % CQdim;
		B_cDim = RQidx2Blk[RQoff]->Cnum;
		E_off = RQidx2Blk[RQoff]->offset + (RQidx2Off[RQoff] * B_cDim) + CQidx2Off[CQoff];
		sB_rDim = RQidx2Dim[RQoff];
		sB_cDim = CQidx2Dim[CQoff];

		int sign01 = 0;
		for(int i = 0; i < swaps.size(); i++)
			sign01 ^= (bonds[swaps[i].b1].Qnums[Q_idxs[swaps[i].b1]].prtF() & bonds[swaps[i].b2].Qnums[Q_idxs[swaps[i].b2]].prtF());
		sign = sign01 ? -1 : 1;

		for(sB_r = 0; sB_r < sB_rDim; sB_r++)
			for(sB_c = 0; sB_c < sB_cDim; sB_c++)
				elem[E_off + (sB_r * B_cDim) + sB_c] *= sign;
	}
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::addGate(std::vector<_Swap>&):");
  }
}

UniTensor& UniTensor::permute(int rowBondNum){
  try{
    std::vector<int> ori_labels = labels;
	  this->permute(ori_labels, rowBondNum);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::permute(int):");
  }
	return *this;
}
UniTensor& UniTensor::permute(int* newLabels, int rowBondNum){
  try{
    std::vector<int> _labels(newLabels, newLabels + bonds.size());
    this->permute(_labels, rowBondNum);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::permute(int*, int):");
  }
	return *this;
}

UniTensor& UniTensor::permute(const std::vector<int>& newLabels, int rowBondNum){
  try{
    if((status & HAVEBOND) == 0){
      std::ostringstream err;
      err<<"There is no bond in the tensor(scalar) to permute.";
      throw std::runtime_error(exception_msg(err.str()));
    }
    if((labels.size() == newLabels.size()) == 0){
      std::ostringstream err;
      err<<"The size of the input new labels does not match for the number of bonds.";
      throw std::runtime_error(exception_msg(err.str()));
    }
    int bondNum = bonds.size();
    std::vector<int> rsp_outin(bondNum);
    int cnt = 0;
    for(int i = 0; i < bondNum; i++)
      for(int j = 0; j < bondNum; j++)
        if(labels[i] == newLabels[j]){
          rsp_outin[j] = i;
          cnt++;
        }
    if((cnt == newLabels.size()) == 0){
      std::ostringstream err;
      err<<"The input new labels do not 1-1 correspond to the labels of the tensor.";
      throw std::runtime_error(exception_msg(err.str()));
    }
    bool inorder = true;
    for(int i = 1; i < bondNum; i++)
      if(((rsp_outin[i] + bondNum - i) % bondNum) != rsp_outin[0]){
        inorder = false;
        break;
      }
    if(inorder && rsp_outin[0] == 0 && RBondNum == rowBondNum)	//do nothing
      return *this;
    else{
      std::vector<Bond> outBonds;
      bool withoutSymmetry = true;
      for(int b = 0; b < bonds.size(); b++){
        outBonds.push_back(bonds[rsp_outin[b]]);
        if(bonds[b].Qnums.size() != 1)
          withoutSymmetry = false;
      }
      for(int b = 0; b < bonds.size(); b++){
        if(b < rowBondNum)
          outBonds[b].change(BD_IN);
        else
          outBonds[b].change(BD_OUT);
      }
      UniTensor UniTout(outBonds, name);
      if(status & HAVEELEM){
        if(withoutSymmetry){
          if(!inorder){
            if(ongpu && UniTout.ongpu){
              size_t* perInfo = (size_t*)malloc(bondNum * 2 * sizeof(size_t));
              std::vector<size_t> newAcc(bondNum);
              newAcc[bondNum - 1] = 1;
              perInfo[bondNum - 1] = 1;
              for(int b = bondNum - 1; b > 0; b--){
                newAcc[b - 1] = newAcc[b] * UniTout.bonds[b].Qdegs[0];
                perInfo[b - 1] = perInfo[b] * bonds[b].Qdegs[0];
              }
              for(int b = 0; b < bondNum; b++)
                perInfo[bondNum + rsp_outin[b]] = newAcc[b];
              double* des_elem = UniTout.elem;
              double* src_elem = elem;
              reshapeElem(src_elem, bondNum, m_elemNum, perInfo, des_elem);
              free(perInfo);
            }
            else{
              double* des_elem = UniTout.elem;
              double* src_elem = elem;
              size_t memsize = m_elemNum * sizeof(DOUBLE);
              if(ongpu){
                src_elem = (double*)elemAllocForce(memsize, false);
                elemCopy(src_elem, elem, memsize, false, ongpu);
              }
              if(UniTout.ongpu)
                des_elem = (double*)elemAllocForce(memsize, false);

              std::vector<size_t> transAcc(bondNum);
              std::vector<size_t> newAcc(bondNum);
              transAcc[bondNum - 1] = 1;
              newAcc[bondNum - 1] = 1;
              for(int b = bondNum - 1; b > 0; b--)
                newAcc[b - 1] = newAcc[b] * UniTout.bonds[b].Qdegs[0];
              std::vector<int> bondDims(bondNum);
              std::vector<int> idxs(bondNum);
              for(int b = 0; b < bondNum; b++){
                transAcc[rsp_outin[b]] = newAcc[b];
                bondDims[b] = bonds[b].Qdegs[0];
                idxs[b] = 0;
              }
              size_t cnt_ot = 0;
              for(int i = 0; i < m_elemNum; i++){
                des_elem[cnt_ot] = src_elem[i];
                for(int bend = bondNum - 1; bend >= 0; bend--){
                  idxs[bend]++;
                  if(idxs[bend] < bondDims[bend]){
                    cnt_ot += transAcc[bend];
                    break;
                  }
                  else{
                    cnt_ot -= transAcc[bend] * (idxs[bend] - 1);
                    idxs[bend] = 0;
                  }
                }
              }
              if(ongpu)
                elemFree(src_elem, memsize, false);
              if(UniTout.ongpu){
                elemCopy(UniTout.elem, des_elem, memsize, UniTout.ongpu, false);
                elemFree(des_elem, memsize, false);
              }
            }
          }
          else{  //non-symmetry inorder
            size_t memsize = m_elemNum * sizeof(DOUBLE);
            elemCopy(UniTout.elem, elem, memsize, UniTout.ongpu, ongpu);
          }
        }
        else{
          int sign = 1;
          //For Fermionic system
          std::vector<_Swap> swaps;
          if(Qnum::isFermionic()){
            std::vector<int> inLabelF(bondNum);
            std::vector<int> outLabelF(bondNum);
            std::vector<int> ordF(bondNum);

            for(int b = 0; b < RBondNum; b++){
              inLabelF[b] = labels[b];
              ordF[b] = b;
            }
            for(int b = 0; b < UniTout.RBondNum; b++)
              outLabelF[b] = newLabels[b];
            for(int b = bondNum - 1; b >= RBondNum; b--){
              ordF[b] = bondNum - b + RBondNum - 1;
              inLabelF[ordF[b]] = labels[b];
            }
            for(int b = bondNum - 1; b >= UniTout.RBondNum; b--)
              outLabelF[bondNum - b + UniTout.RBondNum - 1] = newLabels[b];

            std::vector<int> rspF_outin(bondNum);
            for(int i = 0; i < bondNum; i++)
              for(int j = 0; j < bondNum; j++)
                if(inLabelF[i] == outLabelF[j])
                  rspF_outin[j] = i;
            swaps = recSwap(rspF_outin, ordF);
          }
          //End Fermionic system
          std::vector<int> Qin_idxs(bondNum, 0);
          std::vector<int> Qot_idxs(bondNum, 0);
          int Qin_off, Qot_off;
          int tmp;
          int Qin_RQoff, Qin_CQoff;
          int Qot_CQoff, Qot_RQoff;
          size_t sBin_r, sBin_c;	//sub-block of a Qidx
          size_t sBin_rDim, sBin_cDim;	//sub-block of a Qidx
          size_t sBot_cDim;	//sub-block of a Qidx
          size_t sBot_r, sBot_c;
          size_t Bin_cDim, Bot_cDim;
          size_t Ein_off, Eot_off;
          std::vector<int> sBin_idxs(bondNum, 0);
          std::vector<int> sBin_sBdims(bondNum, 0);
          std::vector<int> Qot_acc(bondNum, 1);
          std::vector<int> sBot_acc(bondNum, 1);
          for(int b = bondNum	- 1; b > 0; b--)
            Qot_acc[b - 1] = Qot_acc[b] * UniTout.bonds[b].Qnums.size();

          for(std::map<int, size_t>::iterator it = QidxEnc.begin(); it != QidxEnc.end(); it++){
            Qin_off = it->first;
            tmp = Qin_off;
            int qdim;
            for(int b = bondNum - 1; b >= 0; b--){
              qdim = bonds[b].Qnums.size();
              Qin_idxs[b] = tmp % qdim;
              sBin_sBdims[b] = bonds[b].Qdegs[Qin_idxs[b]];
              tmp /= qdim;
            }
            Qot_off = 0;
            for(int b = 0; b < bondNum; b++){
              Qot_idxs[b] = Qin_idxs[rsp_outin[b]];
              Qot_off += Qot_idxs[b] * Qot_acc[b];
            }
            for(int b = bondNum	- 1; b > 0; b--)
              sBot_acc[rsp_outin[b-1]] = sBot_acc[rsp_outin[b]] * bonds[rsp_outin[b]].Qdegs[Qot_idxs[b]];
            Qin_RQoff = Qin_off / CQdim;
            Qin_CQoff = Qin_off % CQdim;
            Qot_RQoff = Qot_off / UniTout.CQdim;
            Qot_CQoff = Qot_off % UniTout.CQdim;
            Bin_cDim = RQidx2Blk[Qin_RQoff]->Cnum;
            Bot_cDim = UniTout.RQidx2Blk[Qot_RQoff]->Cnum;
            Ein_off = RQidx2Blk[Qin_RQoff]->offset + (RQidx2Off[Qin_RQoff] * Bin_cDim) + CQidx2Off[Qin_CQoff];
            Eot_off = UniTout.RQidx2Blk[Qot_RQoff]->offset + (UniTout.RQidx2Off[Qot_RQoff] * Bot_cDim) + UniTout.CQidx2Off[Qot_CQoff];
            sBin_rDim = RQidx2Dim[Qin_RQoff];
            sBin_cDim = CQidx2Dim[Qin_CQoff];
            sBot_cDim = UniTout.CQidx2Dim[Qot_CQoff];
            int cnt_ot = 0;
            sBin_idxs.assign(bondNum, 0);
            if(Qnum::isFermionic()){
              int sign01 = 0;
              for(int i = 0; i < swaps.size(); i++)
                sign01 ^= (bonds[swaps[i].b1].Qnums[Qin_idxs[swaps[i].b1]].prtF() & bonds[swaps[i].b2].Qnums[Qin_idxs[swaps[i].b2]].prtF());
              sign = sign01 ? -1 : 1;
            }
            for(sBin_r = 0; sBin_r < sBin_rDim; sBin_r++)
              for(sBin_c = 0; sBin_c < sBin_cDim; sBin_c++){
                sBot_r = cnt_ot / sBot_cDim;
                sBot_c = cnt_ot % sBot_cDim;
                UniTout.elem[Eot_off + (sBot_r * Bot_cDim) + sBot_c] = sign * elem[Ein_off + (sBin_r * Bin_cDim) + sBin_c];
                for(int bend = bondNum - 1; bend >= 0; bend--){
                  sBin_idxs[bend]++;
                  if(sBin_idxs[bend] < sBin_sBdims[bend]){
                    cnt_ot += sBot_acc[bend];
                    break;
                  }
                  else{
                    cnt_ot -= sBot_acc[bend] * (sBin_idxs[bend] - 1);
                    sBin_idxs[bend] = 0;
                  }
                }
              }
          }
        }
        UniTout.status |= HAVEELEM;
      }
      *this = UniTout;
      this->setLabel(newLabels);
    }
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::permute(std::vector<int>&, int):");
  }
  return *this;
}

void UniTensor::setRawElem(const std::vector<DOUBLE>& rawElem){
  try{
    setRawElem(&rawElem[0]);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::setRawElem(std::vector<DOUBLE>&):");
  }
}
void UniTensor::setRawElem(const DOUBLE* rawElem){
  try{
    if((status & HAVEBOND) == 0){
      std::ostringstream err;
      err<<"Setting elements to a tensor without bonds is not supported.";
      throw std::runtime_error(exception_msg(err.str()));
    }
    int bondNum = bonds.size();
    std::vector<int> Q_idxs(bondNum, 0);
    std::vector<int> Q_Bdims(bondNum, 0);
    std::vector<int> sB_idxs(bondNum, 0);
    std::vector<int> sB_sBdims(bondNum, 0);
    std::vector<int> rAcc(bondNum, 1);
    for(int b = 0; b < bondNum; b++)
      Q_Bdims[b] = bonds[b].Qnums.size();
    for(int b = bondNum - 1; b > 0; b--)
      rAcc[b - 1] = rAcc[b] * bonds[b].dim();
    int Q_off;
    int tmp;
    int RQoff, CQoff;
    size_t sB_r, sB_c;	//sub-block of a Qidx
    size_t sB_rDim, sB_cDim;	//sub-block of a Qidx
    size_t B_cDim;
    size_t E_off;
    int R_off;
    double* work = elem;
    if(ongpu){
      work = (double*)malloc(m_elemNum * sizeof(double));
    }
    for(std::map<int, size_t>::iterator it = QidxEnc.begin(); it != QidxEnc.end(); it++){
      Q_off = it->first;
      tmp = Q_off;
      for(int b = bondNum - 1; b >= 0; b--){
        Q_idxs[b] = tmp % Q_Bdims[b];
        tmp /= Q_Bdims[b];
      }
      R_off = 0;
      for(int b = 0; b < bondNum; b++){
        R_off += rAcc[b] * bonds[b].offsets[Q_idxs[b]];
        sB_sBdims[b] = bonds[b].Qdegs[Q_idxs[b]];
      }
      RQoff = Q_off / CQdim;
      CQoff = Q_off % CQdim;
      B_cDim = RQidx2Blk[RQoff]->Cnum;
      E_off = RQidx2Blk[RQoff]->offset + (RQidx2Off[RQoff] * B_cDim) + CQidx2Off[CQoff];
      sB_rDim = RQidx2Dim[RQoff];
      sB_cDim = CQidx2Dim[CQoff];
      sB_idxs.assign(bondNum, 0);
      for(sB_r = 0; sB_r < sB_rDim; sB_r++)
        for(sB_c = 0; sB_c < sB_cDim; sB_c++){
          work[E_off + (sB_r * B_cDim) + sB_c] = rawElem[R_off];
          for(int bend = bondNum - 1; bend >= 0; bend--){
            sB_idxs[bend]++;
            if(sB_idxs[bend] < sB_sBdims[bend]){
              R_off += rAcc[bend];
              break;
            }
            else{
              R_off -= rAcc[bend] * (sB_idxs[bend] - 1);
              sB_idxs[bend] = 0;
            }
          }
        }
    }
    if(ongpu){
      elemCopy(elem, work, m_elemNum * sizeof(double), ongpu, false);
      free(work);
    }
    status |= HAVEELEM;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::setRawElem(double*):");
  }
}

double UniTensor::at(const std::vector<int>& idxs)const{
  try{
    std::vector<size_t> _idxs(idxs.size());
    for(int i = 0; i < idxs.size(); i++)
      _idxs[i] = idxs[i];
    return at(_idxs);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::at(std::vector<int>&):");
    return 0;
  }
}
double UniTensor::at(const std::vector<size_t>& idxs)const{
  try{
    if((status & HAVEBOND) == 0){
      std::ostringstream err;
      err<<"The tensor is a scalar. Use UniTensor::operator[] instead.";
      throw std::runtime_error(exception_msg(err.str()));
    }
    if(!(idxs.size() == bonds.size())){
      std::ostringstream err;
      err<<"The size of input indices array does not match with the number of the bonds.";
      throw std::runtime_error(exception_msg(err.str()));
    }

    int bondNum = bonds.size();
    std::vector<int> Qidxs(bondNum, 0);
    for(int b = 0; b < bondNum; b++){
      if(!(idxs[b] < bonds[b].dim())){
        std::ostringstream err;
        err<<"The input indices are out of range.";
        throw std::runtime_error(exception_msg(err.str()));
      }
      for(int q = bonds[b].offsets.size() - 1; q >= 0; q--){
        if(idxs[b] < bonds[b].offsets[q])
          continue;
        Qidxs[b] = q;
        break;
      }
    }
    std::vector<int> Q_acc(bondNum, 1);
    for(int b = bondNum	- 1; b > 0; b--)
      Q_acc[b - 1] = Q_acc[b] * bonds[b].Qnums.size();
    int Qoff = 0;
    for(int b = 0; b < bondNum; b++)
      Qoff += Q_acc[b] * Qidxs[b];
    if(QidxEnc.find(Qoff) != QidxEnc.end()){
      int Q_RQoff = Qoff / CQdim;
      int Q_CQoff = Qoff % CQdim;
      Block* blk = RQidx2Blk.find(Q_RQoff)->second;
      size_t B_cDim = blk->Cnum;
      size_t sB_cDim = CQidx2Dim.find(Q_CQoff)->second;
      size_t blkRoff = RQidx2Off.find(Q_RQoff)->second;
      size_t blkCoff = CQidx2Off.find(Q_CQoff)->second;
      int boff = blk->offset + (blkRoff * B_cDim) + blkCoff;
      int cnt = 0;
      std::vector<int> D_acc(bondNum, 1);
      for(int b = bondNum	- 1; b > 0; b--)
        D_acc[b - 1] = D_acc[b] * bonds[b].Qdegs[Qidxs[b]];
      for(int b = 0; b < bondNum; b++)
        cnt += (idxs[b] - bonds[b].offsets[Qidxs[b]]) * D_acc[b];
      return elem[boff + (cnt / sB_cDim) * B_cDim + cnt % sB_cDim];
    }
    else{
      return 0.0;;
    }
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::at(std::vector<size_t>&):");
    return 0;
  }
}
double UniTensor::operator[](size_t idx){
  try{
    if(!(idx < m_elemNum)){
      std::ostringstream err;
      err<<"Index exceeds the number of elements("<<m_elemNum<<").";
      throw std::runtime_error(exception_msg(err.str()));
    }
    return getElemAt(idx, elem, ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::operator[](size_t):");
    return 0;
  }
}

UniTensor& UniTensor::transpose(){
  try{
    if(!(status & HAVEBOND)){
      std::ostringstream err;
      err<<"There is no bond in the tensor(scalar) to perform transposition.";
      throw std::runtime_error(exception_msg(err.str()));
    }
    int bondNum = bonds.size();
    std::vector<int> rsp_outin(bondNum);
    int rbondNum = 0;
    for(int b = 0; b < bondNum; b++)
      if(bonds[b].type() == BD_IN)
        rbondNum++;
      else
        break;
    int cbondNum = bondNum - rbondNum;
    for(int b = 0; b < bondNum; b++)
      if(b < cbondNum)
        rsp_outin[b] = rbondNum + b;
      else
        rsp_outin[b] = b - cbondNum;
    std::vector<int> outLabels(bondNum, 0);
    std::vector<Bond> outBonds;
    for(int b = 0; b < bonds.size(); b++){
      outBonds.push_back(bonds[rsp_outin[b]]);
      outLabels[b] = labels[rsp_outin[b]];
    }
    for(int b = 0; b < bondNum; b++){
      if(b < cbondNum)
        outBonds[b].m_type = BD_IN;
      else
        outBonds[b].m_type = BD_OUT;
    }
    UniTensor UniTout(outBonds, name);
    //if(status & HAVELABEL)
    UniTout.setLabel(outLabels);
    if(status & HAVEELEM){
      std::map<Qnum,Block>::iterator it_in;
      std::map<Qnum,Block>::iterator it_out;
      double* elem_in;
      double* elem_out;
      size_t Rnum, Cnum;
      for ( it_in = blocks.begin() ; it_in != blocks.end(); it_in++ ){
        it_out = UniTout.blocks.find((it_in->first));
        Rnum = it_in->second.Rnum;
        Cnum = it_in->second.Cnum;
        elem_in = it_in->second.elem;
        elem_out = it_out->second.elem;

        setTranspose(elem_in, Rnum, Cnum, elem_out, ongpu);
      }
      UniTout.status |= HAVEELEM;
    }
    *this = UniTout;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::transpose():");
  }
  return *this;
}

bool UniTensor::elemCmp(const UniTensor& UniT)const{
  try{
    double diff;
    if(m_elemNum == UniT.m_elemNum){
      for(size_t i = 0; i < m_elemNum; i++){
        diff = fabs(elem[i] - UniT.elem[i]);
        if(diff > 1E-6)
          return false;
      }
    }
    else
      return false;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::elemCmp(uni10::UniTensor&):");
  }
	return true;
}
double UniTensor::trace()const{
  try{
    if(!(status & HAVEELEM)){
      std::ostringstream err;
      err<<"Cannot trace a tensor before setting its elements.";
      throw std::runtime_error(exception_msg(err.str()));
    }
    if(status & HAVEBOND){
      size_t Rnum;
      DOUBLE trVal = 0;
      for(std::map<Qnum, Block>::const_iterator it = blocks.begin() ; it != blocks.end(); it++ ){
        if(!(it->second.Rnum == it->second.Cnum)){
          std::ostringstream err;
          err<<"Cannot trace a non-square block.";
          throw std::runtime_error(exception_msg(err.str()));
        }
        Rnum = it->second.Rnum;
        for(size_t r = 0; r < Rnum; r++)
          trVal += it->second.elem[r * Rnum + r];
      }
      return trVal;
    }
    else
      return elem[0];
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::trace():");
    return 0;
  }
}

UniTensor& UniTensor::partialTrace(int la, int lb){
  try{
    if(!(status & HAVEELEM)){
      std::ostringstream err;
      err<<"Cannot trace bonds of a tensor before setting its elements.";
      throw std::runtime_error(exception_msg(err.str()));
    }
    if(!(bonds.size() > 2)){
      std::ostringstream err;
      err<<"The number of bonds must larger than 2 for performing partialTrace.";
      throw std::runtime_error(exception_msg(err.str()));
    }
    int bondNum = bonds.size();
    std::vector<Bond> newBonds;
    std::vector<int>newLabels(bondNum - 2, 0);
    std::vector<int>rsp_labels(bondNum);
    int ia, ib;
    int enc = 0;
    for(int l = 0; l < labels.size(); l++){
      if(labels[l] == la)
        ia = l;
      else if(labels[l] == lb)
        ib = l;
      else{
        newBonds.push_back(bonds[l]);
        newLabels[enc] = labels[l];
        rsp_labels[enc] = labels[l];
        enc++;
      }
    }
    if(!(enc == newLabels.size())){
      std::ostringstream err;
      err<<"Cannot find the two bonds with the given two labels.";
      throw std::runtime_error(exception_msg(err.str()));
    }

    UniTensor Tt(newBonds, newLabels);
    rsp_labels[bondNum - 2] = labels[ia];
    rsp_labels[bondNum - 1] = labels[ib];
    ia = bondNum - 2;
    ib = bondNum - 1;
    this->permute(rsp_labels, Tt.RBondNum);
    std::vector<int> Q_acc(bondNum, 1);
    for(int b = bondNum - 1; b > 0; b--)
      Q_acc[b - 1] = Q_acc[b] * bonds[b].Qnums.size();
    int tQdim = bonds[ia].Qnums.size();
    /*Sanity Check*/
    if(tQdim == bonds[ib].Qnums.size()){
      std::ostringstream err;
      err<<"The bonds of the given two labels does not match for trace.";
      throw std::runtime_error(exception_msg(err.str()));
    }
    Qnum q0(0, PRT_EVEN);
    for(int q = 0; q < tQdim; q++){
      if(!((bonds[ia].Qnums[q] * bonds[ib].Qnums[q] == q0) && (bonds[ia].Qdegs[q] == bonds[ib].Qdegs[q]))){
        std::ostringstream err;
        err<<"The bonds of the given two labels does not match for trace.";
        throw std::runtime_error(exception_msg(err.str()));
      }
    }
    /*END*/
    int tBnum = Tt.bonds.size();
    //std::vector<int> Q_idxs(bondNum, 0);
    std::vector<int> Qt_Bdims(tBnum, 0);
    for(int b = 0; b < tBnum; b++)
      Qt_Bdims[b] = Tt.bonds[b].Qnums.size();

    int Qt_off;
    int Q_off;
    //int tmp;
    int Qt_RQoff, Qt_CQoff;
    int Q_RQoff, Q_CQoff;
    size_t sBt_rDim, sBt_cDim;	//sub-block of a Qidx of Tt
    size_t sB_rDim, sB_cDim;	//sub-block of a Qidx
    size_t Bt_cDim;
    //int B_cDim;
    size_t Et_off;
    std::vector<size_t> E_offs(tQdim);
    std::vector<size_t> B_cDims(tQdim);
    int tQdim2 = tQdim * tQdim;
    int Qenc = Q_acc[ia] + Q_acc[ib];
    for(std::map<int, size_t>::iterator it = Tt.QidxEnc.begin(); it != Tt.QidxEnc.end(); it++){
      Qt_off = it->first;

      Qt_RQoff = Qt_off / Tt.CQdim;
      Qt_CQoff = Qt_off % Tt.CQdim;
      Bt_cDim = Tt.RQidx2Blk[Qt_RQoff]->Cnum;
      Et_off = Tt.RQidx2Blk[Qt_RQoff]->offset + (Tt.RQidx2Off[Qt_RQoff] * Bt_cDim) + Tt.CQidx2Off[Qt_CQoff];
      sBt_rDim = Tt.RQidx2Dim[Qt_RQoff];
      sBt_cDim = Tt.CQidx2Dim[Qt_CQoff];

      for(int q = 0; q < tQdim; q++){
        Q_off = Qt_off * tQdim2 + q * Qenc;
        Q_RQoff = Q_off / CQdim;
        Q_CQoff = Q_off % CQdim;
        B_cDims[q] = RQidx2Blk[Q_RQoff]->Cnum;
        E_offs[q] = RQidx2Blk[Q_RQoff]->offset + (RQidx2Off[Q_RQoff] * B_cDims[q]) + CQidx2Off[Q_CQoff];
      }
      int tQdeg, sB_c_off;
      DOUBLE trVal;
      for(size_t sB_r = 0; sB_r < sBt_rDim; sB_r++)
        for(size_t sB_c = 0; sB_c < sBt_cDim; sB_c++){
          trVal = 0;
          for(int q = 0; q < tQdim; q++){
            tQdeg = bonds[ia].Qdegs[q];
            sB_c_off = sB_c * (tQdeg * tQdeg);
            for(int t = 0; t < tQdeg; t++){
              trVal += elem[E_offs[q] + (sB_r * B_cDims[q]) + sB_c_off + t * (tQdeg + 1)];
            }
          }
          Tt.elem[Et_off + sB_r * Bt_cDim + sB_c] = trVal;
        }
      Tt.status |= HAVEELEM;
    }
    *this = Tt;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::partialTrace(int, int):");
  }
	return *this;
}
}; /* namespace uni10 */
