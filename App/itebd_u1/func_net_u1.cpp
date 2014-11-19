#include <uni10.hpp>
#include <cstdlib>

//===================================

uni10::UniTensor netLGLGL(uni10::UniTensor ga, uni10::UniTensor gb, uni10::UniTensor la, uni10::UniTensor lb){

        int label_gamma[3];
        int label_lambda[2];
        int label_3b[3];
        int label_4b[4];

        // assign labels to gamma, lambda
        label_lambda[0] = 0;
        label_lambda[1] = 1;
        label_gamma[0] = 1;
        label_gamma[1] = 100;
        label_gamma[2] = 2;

        lb.setLabel(label_lambda);
        ga.setLabel(label_gamma);

        uni10::UniTensor net = uni10::contract(lb, ga, false);
        // false => permutes bonds of original tensors back
        label_3b[0] = 0;
        label_3b[1] = 100;
        label_3b[2] = 2;
        net.permute(label_3b, 1);

        label_lambda[0] = 2;
        label_lambda[1] = 3;
        la.setLabel(label_lambda);
        net = uni10::contract(net, la, false);

        label_3b[0] = 0;
        label_3b[1] = 100;
        label_3b[2] = 3;
        net.permute(label_3b, 1);

	label_gamma[0] = 3;
        label_gamma[1] = 101;
        label_gamma[2] = 4;
        gb.setLabel(label_gamma);

        net = uni10::contract(net, gb, false);
        label_4b[0] = 0;
        label_4b[1] = 100;
        label_4b[2] = 101;
        label_4b[3] = 4;
        net.permute(label_4b, 1);

        label_lambda[0] = 4;
        label_lambda[1] = 5;
        lb.setLabel(label_lambda);
        net = uni10::contract(net, lb, false);

        label_4b[0] = 0;
        label_4b[1] = 100;
        label_4b[2] = 101;
        label_4b[3] = 5;
        net.permute(label_4b, 1);

        return net;
}

//======================================

uni10::UniTensor netLGLGL(uni10::UniTensor gala, uni10::UniTensor gblb, uni10::UniTensor lb){

	int label_3b[3];
        label_3b[0] = gala.label()[2];
        label_3b[1] = gblb.label()[1] + 1;
        label_3b[2] = gblb.label()[2] + 1;
        gblb.setLabel(label_3b);

	uni10::UniTensor net = uni10::contract(gala, gblb, false);

	int label_2b[2];
        label_2b[0] = net.label()[0] - 1;
        label_2b[1] = net.label()[0];
        lb.setLabel(label_2b);
        net = uni10::contract(lb, net, false);

	int label_4b[4] = {0, 100, 101, 5};
        net.setLabel(label_4b);

        return net;
}

//======================================

uni10::UniTensor netLG(uni10::UniTensor ga, uni10::UniTensor la){

        int label_gamma[3];
        int label_lambda[2];
        int label_3b[3];

        // assign labels to gamma, lambda
        label_lambda[0] = 0;
        label_lambda[1] = 1;
        label_gamma[0] = 1;
        label_gamma[1] = 100;
        label_gamma[2] = 2;

        la.setLabel(label_lambda);
        ga.setLabel(label_gamma);

        uni10::UniTensor net = uni10::contract(la, ga, false);
        label_3b[0] = 0;
        label_3b[1] = 100;
        label_3b[2] = 2;
        net.permute(label_3b, 1);

        return net;
}

//=====================================

uni10::UniTensor netGL(uni10::UniTensor ga, uni10::UniTensor la){

        int label_gamma[3];
        int label_lambda[2];
        int label_3b[3];

        // assign labels to gamma, lambda
        label_gamma[0] = 0;
        label_gamma[1] = 100;
        label_gamma[2] = 1;
        label_lambda[0] = 1;
        label_lambda[1] = 2;

        ga.setLabel(label_gamma);
        la.setLabel(label_lambda);

        uni10::UniTensor net = uni10::contract(ga, la, false);
        label_3b[0] = 0;
        label_3b[1] = 100;
        label_3b[2] = 2;
        net.permute(label_3b, 1);

        return net;
}

//======================================

uni10::UniTensor netLGLG(uni10::UniTensor ga, uni10::UniTensor gb, uni10::UniTensor la, uni10::UniTensor lb){

        int label_gamma[3];
        int label_lambda[2];
        int label_3b[3];
	int label_4b[4];

        // assign labels to gamma, lambda
        label_lambda[0] = 0;
        label_lambda[1] = 1;
        label_gamma[0] = 1;
        label_gamma[1] = 100;
        label_gamma[2] = 2;

        lb.setLabel(label_lambda);
        ga.setLabel(label_gamma);
        uni10::UniTensor net = uni10::contract(lb, ga, false);

	label_lambda[0] = 2;
        label_lambda[1] = 3;
	la.setLabel(label_lambda);
	net = uni10::contract(net, la, false);

	label_gamma[0] = 3;
        label_gamma[1] = 101;
        label_gamma[2] = 4;
	gb.setLabel(label_gamma);
	net = uni10::contract(net, gb, false);

	label_4b[0] = 0;
        label_4b[1] = 100;
        label_4b[2] = 101;
	label_4b[3] = 4;
        net.permute(label_4b, 1);

        return net;
}

//================================================

uni10::UniTensor netGLGL(uni10::UniTensor ga, uni10::UniTensor gb, uni10::UniTensor la, uni10::UniTensor lb){

        int label_gamma[3];
        int label_lambda[2];
        int label_3b[3];
        int label_4b[4];

        // assign labels to gamma, lambda
        label_gamma[0] = 0;
        label_gamma[1] = 100;
        label_gamma[2] = 1;
	label_lambda[0] = 1;
        label_lambda[1] = 2;

        ga.setLabel(label_gamma);
	la.setLabel(label_lambda);
        uni10::UniTensor net = uni10::contract(ga, la, false);

	label_gamma[0] = 2;
        label_gamma[1] = 101;
        label_gamma[2] = 3;
        gb.setLabel(label_gamma);
        net = uni10::contract(net, gb, false);

        label_lambda[0] = 3;
        label_lambda[1] = 4;
        lb.setLabel(label_lambda);
        net = uni10::contract(net, lb, false);

        label_4b[0] = 0;
        label_4b[1] = 100;
        label_4b[2] = 101;
        label_4b[3] = 4;
        net.permute(label_4b, 1);

        return net;
}

//===============================================

uni10::UniTensor netGLGL(uni10::UniTensor gala, uni10::UniTensor gblb){

	int label_3b[3];
        label_3b[0] = gala.label()[2];
        label_3b[1] = gblb.label()[1] + 1;
        label_3b[2] = gblb.label()[2] + 1;
        gblb.setLabel(label_3b);

        uni10::UniTensor net = uni10::contract(gala, gblb, false);

        int label_4b[4] = {0, 100, 101, 5};
        net.setLabel(label_4b);
        net.permute(label_4b, 1);

	return net;
}

//===============================================

uni10::UniTensor netMultiCell(int n, uni10::UniTensor ga, uni10::UniTensor gb, uni10::UniTensor la, uni10::UniTensor lb){

	if (n<2){
		std::cout << "number of cells has to be >= 2 for netMultiCell function.\n";
		exit( 1 ); 
	}

	uni10::UniTensor cell = netLGLG(ga, gb ,la, lb);
	uni10::UniTensor net = cell;

	int label_4b[4];

	label_4b[0] = 0;
        label_4b[1] = 100;
        label_4b[2] = 101;
        label_4b[3] = 1;
        net.setLabel(label_4b);

	for (int i=1; i<n-1; i++){

		label_4b[0] = i;
		label_4b[1] = 100 + 2*i;
		label_4b[2] = 101 + 2*i;
		label_4b[3] = i + 1;
		cell.setLabel(label_4b);

		net = uni10::contract(net, cell, false);
	}

	uni10::UniTensor lastCell = netLGLGL(ga, gb, la, lb);
	label_4b[0] = n - 1;
	label_4b[1] = 100 + 2*(n-1);
	label_4b[2] = 101 + 2*(n-1);
	label_4b[3] = n;
	lastCell.setLabel(label_4b);

	net = uni10::contract(net, lastCell, false);
	net.permute(net.label(), 1);

	return net;
}

//===============================================

uni10::UniTensor tranMtx(uni10::UniTensor ket){

	int bn = ket.bondNum();
	int label_nb[bn];

	uni10::UniTensor bra = ket;
	bra.permute(ket.label(), bn-1);

	for (int i=0; i<bn; i++){

		if (i==0 || i==bn-1) label_nb[i] = bra.label()[i] + 1;
		else label_nb[i] = bra.label()[i];
	}

	bra.setLabel(label_nb);

	uni10::UniTensor net = uni10::contract(bra, ket, false);

	int label_4b[4];
	label_4b[0] = net.label()[0];
	label_4b[1] = net.label()[2];
	label_4b[2] = net.label()[1];
	label_4b[3] = net.label()[3];
	net.permute(label_4b, 2);

	return net;
}

//===============================================

/* generate a tensor which is the multiplication of 2^n transfer matrices */
uni10::UniTensor tranMtxChain(uni10::UniTensor ket, int n){

	uni10::UniTensor net = tranMtx(ket);
	uni10::UniTensor cell = tranMtx(ket);

	int label1[4] = {0, 1, 2, 3};
	int label2[4] = {2, 3, 4, 5};

	for (int i=0; i<n; ++i){

		net.setLabel(label1);
		cell.setLabel(label2);

		net = uni10::contract(net, cell, false);
		cell = net;
	}

	return net;
}

//===============================================

uni10::UniTensor tensorOp(uni10::UniTensor ket, uni10::UniTensor op){

	int label_4b[4];

	uni10::UniTensor net;

	if (op.bondNum() == 4){

                int label_op[4];

                label_op[0] = 100;
                label_op[1] = 101;
                label_op[2] = 200;
                label_op[3] = 201;
                op.setLabel(label_op);

                net = uni10::contract(ket, op, false);

                label_4b[0] = 0;
                label_4b[1] = 200;
                label_4b[2] = 201;
                label_4b[3] = 5;
                net.permute(label_4b, 1);
        }
        else if (op.bondNum() == 2){

                int label_op[2];

                label_op[0] = 100;
                label_op[1] = 200;
                op.setLabel(label_op);

                net = uni10::contract(ket, op, false);

                label_4b[0] = 0;
                label_4b[1] = 200;
                label_4b[2] = 101;
                label_4b[3] = 5;
                net.permute(label_4b, 1);
        }

	label_4b[0] = 0;
        label_4b[1] = 100;
        label_4b[2] = 101;
        label_4b[3] = 5;
        net.setLabel(label_4b);

        return net;	
}

//===============================================

uni10::UniTensor expect(uni10::UniTensor ket, uni10::UniTensor op){
	
	int label_4b[4];

	uni10::UniTensor bra = ket;
	bra.permute(ket.label(), ket.bondNum()-1);

	uni10::UniTensor net;

	if (op.bondNum() == 4){

		int label_op[4];

		label_op[0] = 100;
		label_op[1] = 101;
		label_op[2] = 200;
		label_op[3] = 201;
		op.setLabel(label_op);

		net = uni10::contract(ket, op, false);

		label_4b[0] = 0;
		label_4b[1] = 200;
		label_4b[2] = 201;
		label_4b[3] = 5;
		net.permute(label_4b, 1);
	}
	else if (op.bondNum() == 2){

		int label_op[2];

                label_op[0] = 100;
                label_op[1] = 200;
                op.setLabel(label_op);

                net = uni10::contract(ket, op, false);

                label_4b[0] = 0;
                label_4b[1] = 200;
                label_4b[2] = 101;
                label_4b[3] = 5;
                net.permute(label_4b, 1);
	}

	label_4b[0] = 0;
	label_4b[1] = 100;
	label_4b[2] = 101;
	label_4b[3] = 5;
	net.setLabel(label_4b);

	net = uni10::contract(bra, net, false);

	return net;
}

//=============================================

uni10::UniTensor norm(uni10::UniTensor ket){

	uni10::UniTensor bra = ket;
	bra.permute(ket.label(), ket.bondNum()-1);
	uni10::UniTensor net = uni10::contract(bra, ket, false);

	return net;
}

//=============================================


