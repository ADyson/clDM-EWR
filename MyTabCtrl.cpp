// MyTabCtrl.cpp : implementation file
//

#include "stdafx.h"
#include "MyTabCtrl.h"
#include "FTSR.h"
#include "Acquisition.h"
#include "ManualAlignment.h"


// MyTabCtrl

IMPLEMENT_DYNAMIC(MyTabCtrl, CTabCtrl)

//Constructor for the class derived from CTabCtrl
MyTabCtrl::MyTabCtrl()
{
   m_DialogID[0] =IDD_DIALOG1;
   m_DialogID[1] =IDD_DIALOG2;
   m_DialogID[2] = IDD_MANUAL;


   m_Dialog[0] = new EWR();
   m_Dialog[1] = new Acquisition();
   m_Dialog[2] = new ManualAlignment();

   m_nPageCount = 3;

}

//This function creates the Dialog boxes once
void MyTabCtrl::InitDialogs()
{
	m_Dialog[0]->Create(m_DialogID[0],GetParent());
	m_Dialog[1]->Create(m_DialogID[1],GetParent());
	m_Dialog[2]->Create(m_DialogID[2],GetParent());
}



MyTabCtrl::~MyTabCtrl()
{
}


BEGIN_MESSAGE_MAP(MyTabCtrl, CTabCtrl)
	ON_NOTIFY_REFLECT(TCN_SELCHANGE, OnSelchange)
END_MESSAGE_MAP()



// MyTabCtrl message handlers

//Selection change event for the class derived from CTabCtrl
void MyTabCtrl::OnSelchange(NMHDR* pNMHDR, LRESULT* pResult) 
{
   // TODO: Add your control notification handler code here
   ActivateTabDialogs();
   *pResult = 0;
}

void MyTabCtrl::ActivateTabDialogs()
{
  int nSel = GetCurSel();
  if(m_Dialog[nSel]->m_hWnd)
     m_Dialog[nSel]->ShowWindow(SW_HIDE);

  CRect l_rectClient;
  CRect l_rectWnd;

  GetClientRect(l_rectClient);
  AdjustRect(FALSE,l_rectClient);
  GetWindowRect(l_rectWnd);
  GetParent()->ScreenToClient(l_rectWnd);
  l_rectClient.OffsetRect(l_rectWnd.left,l_rectWnd.top);
  for(int nCount=0; nCount < m_nPageCount; nCount++){
     m_Dialog[nCount]->SetWindowPos(&wndTop, l_rectClient.left, l_rectClient.top, l_rectClient.Width(), l_rectClient.Height(), SWP_HIDEWINDOW);
  }
  m_Dialog[nSel]->SetWindowPos(&wndTop, l_rectClient.left, l_rectClient.top, l_rectClient.Width(), l_rectClient.Height(), SWP_SHOWWINDOW);

  m_Dialog[nSel]->ShowWindow(SW_SHOW);

}
