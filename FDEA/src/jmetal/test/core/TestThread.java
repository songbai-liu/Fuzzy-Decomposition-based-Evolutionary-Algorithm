package jmetal.test.core;

import java.io.*;
import java.util.*;

public class TestThread implements Runnable{
	int b = 100;
	
	public synchronized void m1() throws Exception{
		b += 1000;
		Thread.sleep(2000);
		System.out.println(b);
	}
	
	public synchronized void m2() throws Exception{
		Thread.sleep(2000);
		b++;
		System.out.println("m2="+b);
	}
	
	public void run(){
		try{
			m1();

		}catch(Exception e){
			e.printStackTrace();
		}
	}
	
	public static void main(String[] args) throws Exception{
		TestThread tt = new TestThread();
		Thread t1 = new Thread(tt);
		//Thread t2 = new Thread(tt);
		//t1.setPriority(t1.NORM_PRIORITY + 3);
		//t2.setPriority(t2.NORM_PRIORITY - 2);
		t1.start();
		//Thread.sleep(50);
		//t2.start();
		tt.m2();
		Thread.sleep(1);
		Thread.sleep(1);
		Thread.sleep(1);
		Thread.sleep(1);
		Thread.sleep(1);
		System.out.println("tt="+tt.b);
	}

}
