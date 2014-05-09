class ROIPositions
{
public:
	ROIPositions(int T, int L, int B, int R) { iTop = T; iLeft = L; iRight = R; iBottom = B; };
	int iRight;
	int iLeft;
	int iTop;
	int iBottom;

	int Width() { return iRight-iLeft; }
	int Height() { return iBottom-iTop; }
};