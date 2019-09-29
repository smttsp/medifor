#ifndef MAKE_ON_FILTER_H
#define MAKE_ON_FILTER_H

#include <cstring>
#include "Functions.h"

//function f = MakeONFilter(Typef,par)

/* MakeONFilter -- Generate Orthonormal QMF Filter for Wavelet Transform
  Usage
    qmf = MakeONFilter(Typef,par)
  Inputs
    Type   stringf, 'Haar'f, 'Beylkin'f, 'Coiflet'f, 'Daubechies'f,
           'Symmlet'f, 'Vaidyanathan'f,'Battle'
    par    integerf, it is a parameter related to the support and vanishing
           moments of the waveletsf, explained below for each wavelet.

 Outputs
    qmf    quadrature mirror filter

  Description
    The Haar filter (which could be considered a Daubechies-2) was the
    first waveletf, though not called as suchf, and is discontinuous.

    The Beylkin filter places roots for the frequency response function
    close to the Nyquist frequency on the real axis.

    The Coiflet filters are designed to give both the mother and father
    wavelets 2*par vanishing moments; here par may be one of 1f,2f,3f,4 or 5.

    The Daubechies filters are minimal phase filters that generate wavelets
    which have a minimal support for a given number of vanishing moments.
    They are indexed by their lengthf, parf, which may be one of
    4f,6f,8f,10f,12f,14f,16f,18 or 20. The number of vanishing moments is par/2.

    Symmlets are also wavelets within a minimum size support for a given 
    number of vanishing momentsf, but they are as symmetrical as possiblef,
    as opposed to the Daubechies filters which are highly asymmetrical.
    They are indexed by parf, which specifies the number of vanishing
    moments and is equal to half the size of the support. It ranges 
    from 4 to 10.

    The Vaidyanathan filter gives an exact reconstructionf, but does not
    satisfy any moment condition.  The filter has been optimized for
    speech coding.

    The Battle-Lemarie filter generate spline orthogonal wavelet basis.
    The parameter par gives the degree of the spline. The number of 
    vanishing moments is par+1.

  See Also
    FWT_POf, IWT_POf, FWT2_POf, IWT2_POf, WPAnalysis

  References
    The books by Daubechies and Wickerhauser.
	*/


class MakeONFilter_cls {
public:
	float norm(float* arr, int size) {
		float sum = 0;
		for (int i = 0; i < size; i++) {
			sum += arr[i] * arr[i];
		}
		return abs(sqrt(sum));
	}


	Mat MakeONFilter(String type, int par) {
		float tempArr[30]; int size = 0;

		if (type.compare("Haar") == 0) {
			float x = (float)abs(sqrt(0.5));
			float arr[] = { x,x };
			size = 2;
			for (int i = 0; i < size; i++)
				tempArr[i] = arr[i];

		}
		else if (type.compare("Beylkin") == 0) {
			float arr[] = { 0.099305765374f,	0.424215360813f,	0.699825214057f,
							   0.449718251149f, -0.110927598348f, -0.264497231446f,
							   0.026900308804f, 0.155538731877f, -0.017520746267f,
							   -0.088543630623f, 0.019679866044f, 0.042916387274f,
							   -0.017460408696f, -0.014365807969f, 0.010040411845f,
							   0.001484234782f, -0.002736031626f, 0.000640485329f };
			size = sizeof(arr) / sizeof(float);
			for (int i = 0; i < size; i++)
				tempArr[i] = arr[i];
		}

		else if (type.compare("Coiflet") == 0) {
			if (par == 1) {
				float arr[] = { .038580777748f, -0.126969125396f, -0.077161555496f,
								0.607491641386f, 0.745687558934f, 0.226584265197f };
				size = sizeof(arr) / sizeof(float);
				for (int i = 0; i < size; i++)
					tempArr[i] = arr[i];
			}
			else if (par == 2) {
				float arr[] = { 0.016387336463f, -0.041464936782f, -0.067372554722f,
								0.386110066823f, 0.812723635450f, 0.417005184424f,
								-0.076488599078f, -0.059434418646f, 0.023680171947f,
								0.005611434819f, -0.001823208871f, -0.000720549445f };
				size = sizeof(arr) / sizeof(float);
				for (int i = 0; i < size; i++)
					tempArr[i] = arr[i];
			}
			else if (par == 3) {
				float arr[] = { -0.003793512864f, 0.007782596426f, 0.023452696142f,
								-0.065771911281f, -0.061123390003f, 0.405176902410f,
								0.793777222626f, 0.428483476378f, -0.071799821619f,
								-0.082301927106f, 0.034555027573f, 0.015880544864f,
								-0.009007976137f, -0.002574517688f, 0.001117518771f,
								0.000466216960f, -0.000070983303f, -0.000034599773f };
				size = sizeof(arr) / sizeof(float);
				for (int i = 0; i < size; i++)
					tempArr[i] = arr[i];
			}
			else if (par == 4) {
				float arr[] = { 0.000892313668f, -0.001629492013f, -0.007346166328f,
								 0.016068943964f, 0.026682300156f, -0.081266699680f,
								-0.056077313316f, 0.415308407030f, 0.782238930920f,
								0.434386056491f, -0.066627474263f, -0.096220442034f,
								0.039334427123f, 0.025082261845f, -0.015211731527f,
								-0.005658286686f, 0.003751436157f, 0.001266561929f,
								-0.000589020757f, -0.000259974552f, 0.000062339034f,
								0.000031229876f, -0.000003259680f, -0.000001784985f };
				size = sizeof(arr) / sizeof(float);
				for (int i = 0; i < size; i++)
					tempArr[i] = arr[i];
			}
			else if (par == 5) {
				float arr[] = { -0.000212080863f, 0.000358589677f, 0.002178236305f,
								-0.004159358782f, -0.010131117538f, 0.023408156762f,
								0.028168029062f, -0.091920010549f, -0.052043163216f,
								0.421566206729f, 0.774289603740f, 0.437991626228f,
								-0.062035963906f, -0.105574208706f, 0.041289208741f,
								0.032683574283f, -0.019761779012f, -0.009164231153f,
								0.006764185419f, 0.002433373209f, -0.001662863769f,
								-0.000638131296f, 0.000302259520f, 0.000140541149f,
								-0.000041340484f, -0.000021315014f, 0.000003734597f,
								0.000002063806f, -0.000000167408f, -0.000000095158f };
				size = sizeof(arr) / sizeof(float);
				for (int i = 0; i < size; i++)
					tempArr[i] = arr[i];
			}
		}
		else if (type.compare("Daubechies") == 0) {
			if (par == 4) {
				float arr[] = { 0.482962913145f, 0.836516303738f,
								 0.224143868042f, -0.129409522551f };
				size = sizeof(arr) / sizeof(float);
				for (int i = 0; i < size; i++)
					tempArr[i] = arr[i];
			}
			else if (par == 6) {
				float arr[] = { 0.332670552950f, 0.806891509311f,
								0.459877502118f, -0.135011020010f,
								-0.085441273882f, 0.035226291882f };
				size = sizeof(arr) / sizeof(float);
				for (int i = 0; i < size; i++)
					tempArr[i] = arr[i];
			}
			else if (par == 8) {
				float arr[] = { 0.230377813309f, 0.714846570553f,
								 0.630880767930f, -0.027983769417f,
								 -0.187034811719f, 0.030841381836f,
								 0.032883011667f, -0.010597401785f };
				size = sizeof(arr) / sizeof(float);
				for (int i = 0; i < size; i++)
					tempArr[i] = arr[i];
			}
			else if (par == 10) {
				float arr[] = { 0.160102397974f, 0.603829269797f, 0.724308528438f,
								0.138428145901f, -0.242294887066f, -0.032244869585f,
								0.077571493840f, -0.006241490213f, -0.012580751999f,
								0.003335725285f };
				size = sizeof(arr) / sizeof(float);
				for (int i = 0; i < size; i++)
					tempArr[i] = arr[i];
			}
			else if (par == 12) {
				float arr[] = { 0.111540743350f, 0.494623890398f, 0.751133908021f,
								0.315250351709f, -0.226264693965f, -0.129766867567f,
								0.097501605587f, 0.027522865530f, -0.031582039317f,
								0.000553842201f, 0.004777257511f, -0.001077301085f };
				size = sizeof(arr) / sizeof(float);
				for (int i = 0; i < size; i++)
					tempArr[i] = arr[i];
			}
			else if (par == 14) {
				float arr[] = { 0.077852054085f, 0.396539319482f, 0.729132090846f,
								0.469782287405f, -0.143906003929f, -0.224036184994f,
								0.071309219267f, 0.080612609151f, -0.038029936935f,
								-0.016574541631f, 0.012550998556f, 0.000429577973f,
								-0.001801640704f, 0.000353713800f };
				size = sizeof(arr) / sizeof(float);
				for (int i = 0; i < size; i++)
					tempArr[i] = arr[i];
			}
			else if (par == 16) {
				float arr[] = { 0.054415842243f, 0.312871590914f, 0.675630736297f,
								0.585354683654f, -0.015829105256f, -0.284015542962f,
								0.000472484574f, 0.128747426620f, -0.017369301002f,
								-0.044088253931f, 0.013981027917f, 0.008746094047f,
								-0.004870352993f, -0.000391740373f, 0.000675449406f,
								-0.000117476784f };
				size = sizeof(arr) / sizeof(float);
				for (int i = 0; i < size; i++)
					tempArr[i] = arr[i];
			}
			else if (par == 18) {
				float arr[] = { 0.038077947364f, 0.243834674613f, 0.604823123690f,
								0.657288078051f, 0.133197385825f, -0.293273783279f,
								-0.096840783223f, 0.148540749338f, 0.030725681479f,
								-0.067632829061f, 0.000250947115f, 0.022361662124f,
								-0.004723204758f, -0.004281503682f, 0.001847646883f,
								0.000230385764f, -0.000251963189f, 0.000039347320f };
				size = sizeof(arr) / sizeof(float);
				for (int i = 0; i < size; i++)
					tempArr[i] = arr[i];
			}
			else if (par == 20) {
				float arr[] = { 0.026670057901f, 0.188176800078f, 0.527201188932f,
								0.688459039454f, 0.281172343661f, -0.249846424327f,
								-0.195946274377f, 0.127369340336f, 0.093057364604f,
								-0.071394147166f, -0.029457536822f, 0.033212674059f,
								0.003606553567f, -0.010733175483f, 0.001395351747f,
								0.001992405295f, -0.000685856695f, -0.000116466855f,
								0.000093588670f, -0.000013264203f };
				size = sizeof(arr) / sizeof(float);
				for (int i = 0; i < size; i++)
					tempArr[i] = arr[i];
			}

		}
		else if (type.compare("Symmlet") == 0) {
			if (par == 4) {
				float arr[] = { -0.107148901418f, -0.041910965125f, 0.703739068656f,
								1.136658243408f, 0.421234534204f, -0.140317624179f,
								-0.017824701442f, 0.045570345896f };
				size = sizeof(arr) / sizeof(float);
				for (int i = 0; i < size; i++)
					tempArr[i] = arr[i];
			}
			else if (par == 5) {
				float arr[] = { 0.038654795955f, 0.041746864422f, -0.055344186117f,
								0.281990696854f,	1.023052966894f, 0.896581648380f,
								0.023478923136f, -0.247951362613f, -0.029842499869f,
								0.027632152958f };
				size = sizeof(arr) / sizeof(float);
				for (int i = 0; i < size; i++)
					tempArr[i] = arr[i];
			}
			else if (par == 6) {
				float arr[] = { 0.021784700327f, 0.004936612372f, -0.166863215412f,
								-0.068323121587f, 0.694457972958f, 1.113892783926f,
								0.477904371333f, -0.102724969862f, -0.029783751299f,
								0.063250562660f, 0.002499922093f, -0.011031867509f };
				size = sizeof(arr) / sizeof(float);
				for (int i = 0; i < size; i++)
					tempArr[i] = arr[i];
			}
			else if (par == 7) {
				float arr[] = { 0.003792658534f, -0.001481225915f, -0.017870431651f,
								 0.043155452582f, 0.096014767936f, -0.070078291222f,
								 0.024665659489f, 0.758162601964f, 1.085782709814f,
								 0.408183939725f, -0.198056706807f, -0.152463871896f,
								 0.005671342686f, 0.014521394762f };
				size = sizeof(arr) / sizeof(float);
				for (int i = 0; i < size; i++)
					tempArr[i] = arr[i];
			}
			else if (par == 8) {
				float arr[] = { 0.002672793393f, -0.000428394300f, -0.021145686528f,
								0.005386388754f, 0.069490465911f, -0.038493521263f,
								-0.073462508761f, 0.515398670374f, 1.099106630537f,
								0.680745347190f, -0.086653615406f, -0.202648655286f,
								0.010758611751f, 0.044823623042f, -0.000766690896f,
								-0.004783458512f };
				size = sizeof(arr) / sizeof(float);
				for (int i = 0; i < size; i++)
					tempArr[i] = arr[i];
			}
			else if (par == 9) {
				float arr[] = { 0.001512487309f, -0.000669141509f, -0.014515578553f,
								0.012528896242f, 0.087791251554f, -0.025786445930f,
								-0.270893783503f, 0.049882830959f, 0.873048407349f,
								1.015259790832f, 0.337658923602f, -0.077172161097f,
								0.000825140929f, 0.042744433602f, -0.016303351226f,
								-0.018769396836f, 0.000876502539f, 0.001981193736f };
				size = sizeof(arr) / sizeof(float);
				for (int i = 0; i < size; i++)
					tempArr[i] = arr[i];
			}
			else if (par == 10) {
				float arr[] = { 0.001089170447f, 0.000135245020f, -0.012220642630f,
								 -0.002072363923f, 0.064950924579f, 0.016418869426f,
								 -0.225558972234f, -0.100240215031f, 0.667071338154f,
								 1.088251530500f, 0.542813011213f, -0.050256540092f,
								 -0.045240772218f, 0.070703567550f, 0.008152816799f,
								 -0.028786231926f, -0.001137535314f, 0.006495728375f,
								 0.000080661204f, -0.000649589896f };
				size = sizeof(arr) / sizeof(float);
				for (int i = 0; i < size; i++)
					tempArr[i] = arr[i];
			}
		}
		else if (type.compare("Vaidyanathan") == 0) {
			float arr[] = { -0.000062906118f, 0.000343631905f, -0.000453956620f,
				 -0.000944897136f, 0.002843834547f, 0.000708137504f,
				 -0.008839103409f, 0.003153847056f, 0.019687215010f,
				 -0.014853448005f, -0.035470398607f, 0.038742619293f,
				 0.055892523691f, -0.077709750902f, -0.083928884366f,
				 0.131971661417f, 0.135084227129f, -0.194450471766f,
				 -0.263494802488f, 0.201612161775f, 0.635601059872f,
				 0.572797793211f, 0.250184129505f, 0.045799334111f };

			size = sizeof(arr) / sizeof(float);
			for (int i = 0; i < size; i++)
				tempArr[i] = arr[i];
		}
		else if (type.compare("Battle") == 0) {
			if (par == 1) {
				float arr2[] = { 0.578163f, 0.280931f, -0.0488618f, -0.0367309f,
								0.012003f, 0.00706442f, -0.00274588f, -0.00155701f,
								0.000652922f, 0.000361781f, -0.000158601f, -0.0000867523f };
			}

			else if (par == 3) {
				float arr2[] = { 0.541736f, 0.30683f, -0.035498f, -0.0778079f,
								 0.0226846f, 0.0297468f, -0.0121455f,-0.0127154f,
								 0.00614143f,0.00579932f, -0.00307863f,-0.00274529f,
								 0.00154624f, 0.00133086f, -0.000780468f,-0.00065562f,
								 0.000395946f, 0.000326749f, -0.000201818f,-0.000164264f,
								 0.000103307f };
			}

			else if (par == 5) {
				float arr2[] = { 0.528374f, 0.312869f,-0.0261771f, -0.0914068f,
								0.0208414f, 0.0433544f, -0.0148537f, -0.0229951f,
								0.00990635f, 0.0128754f, -0.00639886f, -0.00746848f,
								0.00407882f, 0.00444002f, -0.00258816f, -0.00268646f,
								0.00164132f, 0.00164659f, -0.00104207f, -0.00101912f,
								0.000662836f, 0.000635563f, -0.000422485f, -0.000398759f,
								0.000269842f, 0.000251419f, -0.000172685f, -0.000159168f,
								0.000110709f, 0.000101113f };
			}
		}
		Mat qmf;
		qmf.create(Size(size, 1), CV_32FC1);
		float normRes = norm(tempArr, size);

		for (int i = 0; i < size; i++) {
			qmf.at<float>(0, i) = tempArr[i] / normRes;
		}
		return qmf;
	}
};
#endif