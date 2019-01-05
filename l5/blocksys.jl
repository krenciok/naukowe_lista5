module blocksys
#=
Tomasz Krent 236595
=#
	export readA,readB,writeB,z1a,z1b,z2a,z2b,z3,matcond,blockmat,zadanie1a,zadanie1b,zadanie2a,zadanie2b
	
	#READA
	#
	#
	function readA(inputfile::String)
		m=Array{Float64}
		m=readdlm(inputfile)
		n=Float64(m[1,1])
		l=Float64(m[1,2])
		v=n/l
		liczba=l*l*v + l*(v-1)*2
		
		Wx=Array{Float64}(Int(n),Int(l)*2+2)
		Wx=zeros(Int(n),Int(l)*2+2)
		
		
		for i in 2:1:Int(liczba+1)
			a=m[i,1]
			a1=m[i,2]
			if a<=l
				is=0
			else
				is=1
			end
			if(floor((a-1)/l)==floor((a1-1)/l))
				c=a1%l
				if c==0
					c=l
				end
				Wx[Int(a),Int(c)+is]=m[i,3]
			elseif a<a1
				c=a1%l
				if c==0
					c=l
				end
				Wx[Int(a),Int(l)+is+Int(c)]=m[i,3]
			else
				Wx[Int(a),1]=m[i,3]
			end
		
		end
		return (Wx,Int(n),Int(l))
	end
	#
	#
	#
	function readB(inputfile::String)
		m2=Array{Float64}
		m2=readdlm(inputfile)
		n2=Float64(m2[1,1])
		
		B=Array{Float64}(Int(n2),1)
		
		for i in 2:1:Int(n2)+1
			B[Int(i)-1,1]=m2[i,1]
			
		end
		return B
	end

	#WRITEB
	#
	#
	function writeB(Wx::Array,n::Int,l::Int,outputfile::String)
		B=Array{Float64}(Int(n)+1,1)
		B[1,1]=Int(n)
		for i in 1:1:Int(n)
			B[i+1]=0.0
			mini=min(Int(n)-i+1,Int(l)*2)
			for k in 1:1:Int(l)*2+1
				B[i+1,1]=B[i+1,1]+Wx[i,k]
			end	
		end
		writedlm(outputfile,B)
	end
	
	#
	#
	#
	function z1a(Wx::Array,B::Array,n::Int,l::Int,outputfile::String)

		X=Array{Float64}(Int(n),1)
		
		zm=Int(l)-1
		
		for i in 1:1:Int(n)-1
			for j in 1:1:zm
				tmp=Wx[i+j,1]/Wx[i,1]
				mini=min(Int(n)-i+1,Int(l)*2)
				for k in 1:1:Int(l)*2+1
					Wx[i+j,k]=Wx[i+j,k+1]-tmp*Wx[i,k+1]
				end
				B[i+j,1]= B[i+j,1] - tmp*B[i]
			end
			zm=zm-1
			if zm==0
				zm=Int(l)
			end
		end
		
		k=0
		bl=Array{Float64}(1,1)
		bl[1,1]=0.0
		for i in Int(n):-1:1
			X[i,1]=B[i,1]
			for j in 1:1:k
				X[i,1]=X[i,1]-Wx[i,1+j]*X[i+j,1]
			end
			X[i,1]=X[i,1]/Wx[i,1]
			if k<2*l+1
				k=k+1
			end
			bl[1,1]=bl[1,1]+abs(1.0-X[i,1])
		end
		
		bl[1,1]=bl[1,1]/n
		
		writedlm(outputfile,[bl;X])
		
		return bl[1,1]
		
	end
	
	#
	#
	#
	function z1b(Wx::Array,B::Array,n::Int,l::Int,outputfile::String)

		X=Array{Float64}(Int(n),1)

		zm=Int(l)-1
		
		for i in 1:1:Int(n)-1
			maxi=0
			t=i
			for j in 0:1:zm
				if abs(Wx[i+j,1])>maxi
					maxi=abs(Wx[i+j,1])
					t=i+j
				end
			end
			
			Wx[i,:],Wx[t,:]=Wx[t,:],Wx[i,:]
			B[i,:],B[t,:]=B[t,:],B[i,:]
			
			for j in 1:1:zm
				tmp=Wx[i+j,1]/Wx[i,1]
				mini=min(Int(n)-i+1,Int(l)*2)
				for k in 1:1:Int(l)*2+1
					Wx[i+j,k]=Wx[i+j,k+1]-tmp*Wx[i,k+1]
				end
				B[i+j,1]= B[i+j,1] - tmp*B[i]
			end
			zm=zm-1
			if zm==0
				zm=Int(l)
			end
		end
		
		k=0
		bl=Array{Float64}(1,1)
		bl[1,1]=0.0
		for i in Int(n):-1:1
			X[i,1]=B[i,1]
			for j in 1:1:k
				X[i,1]=X[i,1]-Wx[i,1+j]*X[i+j,1]
			end
			X[i,1]=X[i,1]/Wx[i,1]
			if k<2*l+1
				k=k+1
			end
			bl[1,1]=bl[1,1]+abs(1.0-X[i,1])
		end
		
		bl[1,1]=bl[1,1]/n
		
		writedlm(outputfile,[bl;X])
		
		return bl[1,1]
		
	end
	
	#
	#
	#
	function z2a(Wx::Array,B::Array,n::Int,l::Int)

		Wx2=Array{Float64}(Int(n),Int(l)+1)
		Wx2=zeros(Int(n),Int(l)+1)

		zm=Int(l)-1
		ii=1
		Wx2[1,1]=1.0
		for i in 1:1:Int(n)-1
			for j in 1:1:zm
				tmp=Wx[i+j,1]/Wx[i,1]
				mini=min(Int(n)-i+1,Int(l)*2)
				for k in 1:1:Int(l)*2+1
					Wx[i+j,k]=Wx[i+j,k+1]-tmp*Wx[i,k+1]
					Wx2[i+j,ii]=tmp
				end
				
			end
			
			ii=ii+1
			Wx2[i+1,ii]=1
			zm=zm-1
			if zm==0
				zm=Int(l)
				ii=1
			end	
		end
		return (Wx,Wx2,B)
	end
	
	#
	#
	#
	function z2b(Wx::Array,B::Array,n::Int,l::Int)

		Wx2=Array{Float64}(Int(n),Int(l)+1)
		Wx2=zeros(Int(n),Int(l)+1)

		zm=Int(l)-1
		ii=1
		Wx2[1,1]=1.0
		for i in 1:1:Int(n)-1
		
			maxi=0
			t=i
			for j in 0:1:zm
				if abs(Wx[i+j,1])>maxi
					maxi=abs(Wx[i+j,1])
					t=i+j
				end
			end
			
			Wx[i,:],Wx[t,:]=Wx[t,:],Wx[i,:]
			Wx2[i,:],Wx2[t,:]=Wx2[t,:],Wx2[i,:]
			B[i,:],B[t,:]=B[t,:],B[i,:]
			
			for j in 1:1:zm
				tmp=Wx[i+j,1]/Wx[i,1]
				mini=min(Int(n)-i+1,Int(l)*2)
				for k in 1:1:Int(l)*2+1
					Wx[i+j,k]=Wx[i+j,k+1]-tmp*Wx[i,k+1]
					Wx2[i+j,ii]=tmp
				end
				
			end
			
			ii=ii+1
			Wx2[i+1,ii]=1
			zm=zm-1
			if zm==0
				zm=Int(l)
				ii=1
			end	
		end
		return (Wx,Wx2,B)
		
	end
	
	#
	#
	#
	function z3(Wx::Array,Wx2::Array,B::Array,n::Int,l::Int,outputfile::String)
		X=Array{Float64}(Int(n),1)
		Y=Array{Float64}(Int(n),1)
		
		k=0
		ii=1
		for i in 1:1:Int(n)
			Y[i,1]=B[i,1]
			for j in 1:1:k
				Y[i,1]=Y[i,1]-Wx2[i,k+1-j]*Y[i-j,1]
			end
			k=k+1
			if k+ii==Int(l)+1
				k=1
				ii=0
			end
			
		end
		
		k=0
		bl=Array{Float64}(1,1)
		bl[1,1]=0.0
		for i in Int(n):-1:1
			X[i,1]=Y[i,1]
			for j in 1:1:k
				X[i,1]=X[i,1]-Wx[i,1+j]*X[i+j,1]
			end
			X[i,1]=X[i,1]/Wx[i,1]
			if k<2*l+1
				k=k+1
			end
			bl[1,1]=bl[1,1]+abs(1.0-X[i,1])
		end
		
		bl[1,1]=bl[1,1]/n
		
		writedlm(outputfile,[bl;X])
		
		return bl[1,1]
		
	end

	# Pawel Zielinski

	function matcond(n::Int, c::Float64)
		# Function generates a random square matrix A of size n with
		# a given condition number c.	
		# Inputs:
		#	n: size of matrix A, n>1
		#	c: condition of matrix A, c>= 1.0	
		#
		# Usage: matcond (10, 100.0);
		#
		
        if n < 2
         error("size n should be > 1")
        end
        if c< 1.0
         error("condition number  c of a matrix  should be >= 1.0")
        end
        (U,S,V)=svd(rand(n,n))
        return U*diagm(linspace(1.0,c,n))*V'
	end # matcond

	# Pawel Zielinski
	
	function blockmat(n::Int, l::Int, ck::Float64, outputfile::String)
		# Function generates a random block sparse matrix A of size n with
		# a given condition number ck of inner block Ak and it save the output
		# matrix in a text file.
		# Inputs:
		#	n: size of block matrix A, n>1
		# l: size of inner matrices Ak, n mod l =0 (n is  divisible by l)
		#	ck: condition of inner matrix Ak, ck>= 1.0	
		# outputfile: name of the output text file
		#
		# Usage: blockmat(100, 4 ,10.0, "A.txt")
		#		
		#
		#  the output file format
	  #  n  l              <--- the size of block matrix A, the size of inner matrices Ak
		#  i1  j1   A[i1,j1] <--- a non-zero element of block matrix A 
		#  i2  j2   A[i2,j2] <--- a non-zero element of block matrix A 
		#  i3  j3   A[i3,j3] <--- a non-zero element of block matrix A 
		#  ...
		#  ...
		#  EOF
		#
					
		if n < 2
		 error("size n should be > 1")
		end
		if n%l!=0 
			error("n is not divisible by l")
		end
					
		nb=div(n,l)
		Ak=Array{Float64}(l,l)		
		open(outputfile, "w") do f
			println(f, n," ",l)
			for k in 1:nb
				Ak=matcond(l, ck)
				for i in 1:l, j in 1:l
					println(f,(k-1)*l+i," ",(k-1)*l+j," ", Ak[i,j])
				end
				if k<nb
			  	 for i in 1:l
						println(f,(k-1)*l+i," ",k*l+i," ",0.3*rand())
				 	 end
				end
				if k>1
			   for i in 1:l
					 println(f,(k-1)*l+i," ",(k-1)*l," ",0.3*rand())
				 end
				end 
			end
		end	 # do
	end # blockmat
		
	function zadanie1a()
	
		T1=Array{Float64}(2,49)
		T2=Array{Float64}(2,49)
		T3=Array{Float64}(2,49)
		licz=1
		for it in 100:100:4900
		
			numer=5
			blockmat(it,numer,100.0,"zadanie1a/$(it)_$(numer)_A.txt")
			q1=readA("zadanie1a/$(it)_$(numer)_A.txt")
			writeB(q1[1],q1[2],q1[3],"zadanie1a/$(it)_$(numer)_b.txt")
			q2=readB("zadanie1a/$(it)_$(numer)_b.txt")
			T1[2,licz]=z1a(q1[1],q2,q1[2],q1[3],"zadanie1a/$(it)_$(numer)_x.txt")
			T1[1,licz]=it
			
			numer=10
			blockmat(it,numer,100.0,"zadanie1a/$(it)_$(numer)_A.txt")
			q1=readA("zadanie1a/$(it)_$(numer)_A.txt")
			writeB(q1[1],q1[2],q1[3],"zadanie1a/$(it)_$(numer)_b.txt")
			q2=readB("zadanie1a/$(it)_$(numer)_b.txt")
			T2[2,licz]=z1a(q1[1],q2,q1[2],q1[3],"zadanie1a/$(it)_$(numer)_x.txt")
			T2[1,licz]=it
			
			numer=20
			blockmat(it,numer,100.0,"zadanie1a/$(it)_$(numer)_A.txt")
			q1=readA("zadanie1a/$(it)_$(numer)_A.txt")
			writeB(q1[1],q1[2],q1[3],"zadanie1a/$(it)_$(numer)_b.txt")
			q2=readB("zadanie1a/$(it)_$(numer)_b.txt")
			T3[2,licz]=z1a(q1[1],q2,q1[2],q1[3],"zadanie1a/$(it)_$(numer)_x.txt")
			T3[1,licz]=it
			
			licz=licz+1
		end
		return(T1,T2,T3)
	end
	
	function zadanie1b()
	
		T1=Array{Float64}(2,49)
		T2=Array{Float64}(2,49)
		T3=Array{Float64}(2,49)
		licz=1
		for it in 100:100:4900
		
			numer=5
			blockmat(it,numer,100.0,"zadanie1b/$(it)_$(numer)_A.txt")
			q1=readA("zadanie1b/$(it)_$(numer)_A.txt")
			writeB(q1[1],q1[2],q1[3],"zadanie1b/$(it)_$(numer)_b.txt")
			q2=readB("zadanie1b/$(it)_$(numer)_b.txt")
			T1[2,licz]=z1b(q1[1],q2,q1[2],q1[3],"zadanie1b/$(it)_$(numer)_x.txt")
			T1[1,licz]=it
			
			numer=10
			blockmat(it,numer,100.0,"zadanie1b/$(it)_$(numer)_A.txt")
			q1=readA("zadanie1b/$(it)_$(numer)_A.txt")
			writeB(q1[1],q1[2],q1[3],"zadanie1b/$(it)_$(numer)_b.txt")
			q2=readB("zadanie1b/$(it)_$(numer)_b.txt")
			T2[2,licz]=z1b(q1[1],q2,q1[2],q1[3],"zadanie1b/$(it)_$(numer)_x.txt")
			T2[1,licz]=it
			
			numer=20
			blockmat(it,numer,100.0,"zadanie1b/$(it)_$(numer)_A.txt")
			q1=readA("zadanie1b/$(it)_$(numer)_A.txt")
			writeB(q1[1],q1[2],q1[3],"zadanie1b/$(it)_$(numer)_b.txt")
			q2=readB("zadanie1b/$(it)_$(numer)_b.txt")
			T3[2,licz]=z1b(q1[1],q2,q1[2],q1[3],"zadanie1b/$(it)_$(numer)_x.txt")
			T3[1,licz]=it
			
			licz=licz+1
		end
		return(T1,T2,T3)
	end
 
	function zadanie2a()
	
		T1=Array{Float64}(2,49)
		T2=Array{Float64}(2,49)
		T3=Array{Float64}(2,49)
		licz=1
		for it in 100:100:4900
		
			numer=5
			blockmat(it,numer,100.0,"zadanie2a/$(it)_$(numer)_A.txt")
			q1=readA("zadanie2a/$(it)_$(numer)_A.txt")
			writeB(q1[1],q1[2],q1[3],"zadanie2a/$(it)_$(numer)_b.txt")
			q2=readB("zadanie2a/$(it)_$(numer)_b.txt")
			q3=z2a(q1[1],q2,q1[2],q1[3])
			T1[2,licz]=z3(q3[1],q3[2],q3[3],q1[2],q1[3],"zadanie2a/$(it)_$(numer)_x.txt")
			T1[1,licz]=it
			
			numer=10
			blockmat(it,numer,100.0,"zadanie2a/$(it)_$(numer)_A.txt")
			q1=readA("zadanie2a/$(it)_$(numer)_A.txt")
			writeB(q1[1],q1[2],q1[3],"zadanie2a/$(it)_$(numer)_b.txt")
			q2=readB("zadanie2a/$(it)_$(numer)_b.txt")
			q3=z2a(q1[1],q2,q1[2],q1[3])
			T2[2,licz]=z3(q3[1],q3[2],q3[3],q1[2],q1[3],"zadanie2a/$(it)_$(numer)_x.txt")
			T2[1,licz]=it
			
			numer=20
			blockmat(it,numer,100.0,"zadanie2a/$(it)_$(numer)_A.txt")
			q1=readA("zadanie2a/$(it)_$(numer)_A.txt")
			writeB(q1[1],q1[2],q1[3],"zadanie2a/$(it)_$(numer)_b.txt")
			q2=readB("zadanie2a/$(it)_$(numer)_b.txt")
			q3=z2a(q1[1],q2,q1[2],q1[3])
			T3[2,licz]=z3(q3[1],q3[2],q3[3],q1[2],q1[3],"zadanie2a/$(it)_$(numer)_x.txt")
			T3[1,licz]=it
			
			licz=licz+1
		end
		return(T1,T2,T3)
	end

	function zadanie2b()
	
		T1=Array{Float64}(2,49)
		T2=Array{Float64}(2,49)
		T3=Array{Float64}(2,49)
		licz=1
		for it in 100:100:4900
		
			numer=5
			blockmat(it,numer,100.0,"zadanie2b/$(it)_$(numer)_A.txt")
			q1=readA("zadanie2b/$(it)_$(numer)_A.txt")
			writeB(q1[1],q1[2],q1[3],"zadanie2b/$(it)_$(numer)_b.txt")
			q2=readB("zadanie2b/$(it)_$(numer)_b.txt")
			q3=z2b(q1[1],q2,q1[2],q1[3])
			T1[2,licz]=z3(q3[1],q3[2],q3[3],q1[2],q1[3],"zadanie2b/$(it)_$(numer)_x.txt")
			T1[1,licz]=it
			
			numer=10
			blockmat(it,numer,100.0,"zadanie2b/$(it)_$(numer)_A.txt")
			q1=readA("zadanie2b/$(it)_$(numer)_A.txt")
			writeB(q1[1],q1[2],q1[3],"zadanie2b/$(it)_$(numer)_b.txt")
			q2=readB("zadanie2b/$(it)_$(numer)_b.txt")
			q3=z2b(q1[1],q2,q1[2],q1[3])
			T2[2,licz]=z3(q3[1],q3[2],q3[3],q1[2],q1[3],"zadanie2b/$(it)_$(numer)_x.txt")
			T2[1,licz]=it
			
			numer=20
			blockmat(it,numer,100.0,"zadanie2b/$(it)_$(numer)_A.txt")
			q1=readA("zadanie2b/$(it)_$(numer)_A.txt")
			writeB(q1[1],q1[2],q1[3],"zadanie2b/$(it)_$(numer)_b.txt")
			q2=readB("zadanie2b/$(it)_$(numer)_b.txt")
			q3=z2b(q1[1],q2,q1[2],q1[3])
			T3[2,licz]=z3(q3[1],q3[2],q3[3],q1[2],q1[3],"zadanie2b/$(it)_$(numer)_x.txt")
			T3[1,licz]=it
			
			licz=licz+1
		end
		return(T1,T2,T3)
	end
	
end
