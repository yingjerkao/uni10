uni10::UniTensor netLGLGL(uni10::UniTensor ga, uni10::UniTensor gb, uni10::UniTensor la, uni10::UniTensor lb);
uni10::UniTensor netLGLGL(uni10::UniTensor gala, uni10::UniTensor gblb, uni10::UniTensor lb);
uni10::UniTensor netLG(uni10::UniTensor ga, uni10::UniTensor la);
uni10::UniTensor netGL(uni10::UniTensor ga, uni10::UniTensor la);
uni10::UniTensor netLGLG(uni10::UniTensor ga, uni10::UniTensor gb, uni10::UniTensor la, uni10::UniTensor lb);
uni10::UniTensor netGLGL(uni10::UniTensor ga, uni10::UniTensor gb, uni10::UniTensor la, uni10::UniTensor lb);
uni10::UniTensor netGLGL(uni10::UniTensor gala, uni10::UniTensor gblb);

uni10::UniTensor netMultiCell(int n, uni10::UniTensor ga, uni10::UniTensor gb, uni10::UniTensor la, uni10::UniTensor lb);

uni10::UniTensor tranMtx(uni10::UniTensor ket);
uni10::UniTensor tranMtxChain(uni10::UniTensor ket, int n);

uni10::UniTensor tensorOp(uni10::UniTensor ket, uni10::UniTensor op);
uni10::UniTensor expect(uni10::UniTensor ket, uni10::UniTensor op);
uni10::UniTensor norm(uni10::UniTensor ket);


