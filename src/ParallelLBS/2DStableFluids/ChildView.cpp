
// ChildView.cpp : implementation of the CChildView class
//

#include "stdafx.h"
#include "2DStableFluids.h"
#include "ChildView.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif


// CChildView

CChildView::CChildView()
{
	lattice = new Lattice(30, 30, 2);
	windowSize = 600;
	dx = windowSize/lattice->GetXDIM();
	showGrid = false;
	m_timer = 0;
	leftButton = false;
	rightButton = false;

	frameRate = 0;
	timeSinceFrameRateUpdate = 0;
	framesSinceFrameRateUpdate = 0;

}

CChildView::~CChildView()
{
	delete lattice;
}


BEGIN_MESSAGE_MAP(CChildView, CWnd)
	ON_WM_PAINT()
	ON_WM_TIMER()
	ON_WM_LBUTTONDOWN()
	ON_WM_LBUTTONUP()
	ON_WM_RBUTTONDOWN()
	ON_WM_RBUTTONUP()
	ON_WM_MOUSEMOVE()
	ON_WM_KEYDOWN()
END_MESSAGE_MAP()



// CChildView message handlers

BOOL CChildView::PreCreateWindow(CREATESTRUCT& cs) 
{
	if (!CWnd::PreCreateWindow(cs))
		return FALSE;

	cs.dwExStyle |= WS_EX_CLIENTEDGE;
	cs.style &= ~WS_BORDER;
	cs.lpszClass = AfxRegisterWndClass(CS_HREDRAW|CS_VREDRAW|CS_DBLCLKS, 
		::LoadCursor(NULL, IDC_ARROW), reinterpret_cast<HBRUSH>(COLOR_WINDOW+1), NULL);

	return TRUE;
}

void CChildView::OnPaint() 
{	
	CPaintDC dc(this); // device context for painting
	this->GetClientRect(rect);
    // double buffering to remove flickering
    MemDC.CreateCompatibleDC(NULL);
    MemBitmap.CreateCompatibleBitmap(&dc,windowSize+1,windowSize+1);
    CBitmap *pOldBit=MemDC.SelectObject(&MemBitmap);
	MemDC.FillSolidRect(0, 0, windowSize+1, windowSize+1, RGB(255,255,255));
    // clear the screen
	CBrush brush;
    // repaint the window with solid black background
    brush.CreateSolidBrush(RGB(255,255,255));

	//CPen qCirclePen(PS_SOLID, 5, RGB(0, 0, 0));
    //CPen* pqOrigPen = MemDC.SelectObject(&qCirclePen);

	int gridSize = lattice->GetXDIM();
	int c;
	for (int cell_i = 0; cell_i < gridSize; cell_i++)
		for (int cell_j = 0; cell_j < gridSize; cell_j++)
		{
			c = (int) ((36 - lattice->density(cell_i + cell_j * gridSize)) * 9);
			if (c < 0)
				c = 0;
			if (c > 255)
				c = 255;                		
			MemDC.FillSolidRect(cell_i * dx, cell_j * dx, dx, dx, RGB(c, c, c));
		}

	//Draw Grid lines
	if (showGrid)
	{
		// Primal edges
		CPen qLinePen(PS_SOLID, 1, RGB(125, 125, 125));
		MemDC.SelectObject(&qLinePen);
		for (int cell = 0; cell < gridSize; cell++)
		{
			// hrizental grid lines
			MemDC.MoveTo(0, (cell + 1) * dx);
			MemDC.LineTo(windowSize, (cell + 1) * dx);
			//vertical grid lines
			MemDC.MoveTo((cell + 1) * dx, 0);
			MemDC.LineTo((cell + 1) * dx, windowSize);
		}
	}

	timeSinceFrameRateUpdate += lattice->LastTimeStep();
	framesSinceFrameRateUpdate++;

	if (timeSinceFrameRateUpdate > 1000)
	{
		frameRate = (framesSinceFrameRateUpdate * 1000) / timeSinceFrameRateUpdate;
		framesSinceFrameRateUpdate = 0;
		timeSinceFrameRateUpdate  = 0;
	}

	dc.BitBlt(0,0,windowSize+1,windowSize+1,&MemDC,0,0,SRCCOPY);
    MemBitmap.DeleteObject();
    MemDC.DeleteDC();


	int TextWidth = 250;
	int TextHeight = 300;
    MemDC1.CreateCompatibleDC(NULL);
    MemBitmap1.CreateCompatibleBitmap(&dc,TextWidth,TextHeight);
    CBitmap *pOldBit1=MemDC1.SelectObject(&MemBitmap1);
	//Background color for text
	MemDC1.FillSolidRect(windowSize+1, 0, TextWidth, TextHeight, RGB(0,20,20));
  
	MemDC1.SetTextColor(RGB(0,255,255));
	CString s1, s2;
	s1 = "n = ";
	s2.Format(_T("%d"),gridSize);
	s2 = s1+s2;
	MemDC1.TextOutW(3,10,s2);

	MemDC1.SetTextColor(RGB(255,255,255));
	int row = 35;
	MemDC1.TextOutW(3,row,_T("User Interface Guide:"));
	row += 20;
	MemDC1.TextOutW(8,row,_T("Z : Start Animation"));
	row += 20;
	MemDC1.TextOutW(8,row,_T("Left button: Pour Water"));
	row += 20;
	MemDC1.TextOutW(8,row,_T("R : Reset"));
	row += 20;
	MemDC1.TextOutW(8,row,_T("G : Show/Hide Gridlines"));
	oss << "Frame Rate: " << frameRate;
	frameRateStr = oss.str();
	oss.str(std::string());
	oss.flush();
	row += 60;
	MemDC1.TextOutW(8,row, CString(frameRateStr.c_str()));


	dc.BitBlt(windowSize+1,0,TextWidth,TextHeight,&MemDC1,0,0,NOTSRCCOPY);
    MemBitmap1.DeleteObject();
    MemDC1.DeleteDC();
}

void CChildView::OnTimer(UINT_PTR nIDEvent)
{
	if (leftButton) {
		int index = Find_Cell_Index(current_point);	 
		// Inject density
		lattice->SetDensity(index, 5);
	}

	lattice->update();

	Invalidate(false);
	CWnd::OnTimer(nIDEvent);
}


void CChildView::OnLButtonDown(UINT nFlags, CPoint point)
{
	leftButton = true;
	current_point = old_point = point;
	CWnd::OnLButtonDown(nFlags, point);
}


void CChildView::OnLButtonUp(UINT nFlags, CPoint point)
{
	leftButton = false;
	CWnd::OnLButtonUp(nFlags, point);
}


void CChildView::OnRButtonDown(UINT nFlags, CPoint point)
{
	rightButton = true;
	current_point = old_point = point;
	CWnd::OnRButtonDown(nFlags, point);
}


void CChildView::OnRButtonUp(UINT nFlags, CPoint point)
{
	rightButton = false;

	CWnd::OnRButtonUp(nFlags, point);
}


void CChildView::OnMouseMove(UINT nFlags, CPoint point)
{
	if (leftButton || rightButton) {
		old_point = current_point;
		current_point = point;
	}

	CWnd::OnMouseMove(nFlags, point);
}


void CChildView::OnKeyDown(UINT nChar, UINT nRepCnt, UINT nFlags)
{

	switch( nChar)
	{
	case 'a':
	case 'A':
		lattice->update();
		break;
	case 'Z':
	case 'z':
		if (m_timer) {
			KillTimer(m_timer);
			m_timer = 0;
		}
		else
		{
		   m_timer = SetTimer(1, (int) (25), NULL);
		}
		lattice->SetLastUpdateTime();
		break;
	case 'r':
	case 'R':
		lattice->Reset();
		Invalidate(false);
		break;
	case 'G': // Show/Hide the grid
	case 'g':
		showGrid = !showGrid;
		InvalidateRect(NULL,FALSE);
		break;
	}

	CWnd::OnKeyDown(nChar, nRepCnt, nFlags);
}

int CChildView::Find_Cell_Index(CPoint point)
{
	int x = point.x;
	int y = point.y;
	// set boundaries for mouse input to be inside the rectangular
	if (x >= windowSize)
		x = windowSize - 1;
	if (x < 0)
		x = 0;
	if (y >= windowSize)
		y = windowSize - 1;
	if (y < 0)
		y = 0;
	// finwindowSize the cell inwindowSizeex of mouse input
	int cell_i = x / dx;
	int cell_j = y / dx;

	return cell_i + lattice->GetXDIM() * cell_j;
}
