for test in public-instances/*in; do
   echo "===========> Running test file $test";
   ./main $test > ${test%in}tst;
   diff ${test%in}tst ${test%in}out
   if diff ${test%in}tst ${test%in}out >/dev/null ; then
      echo PASSED!
   else
      echo FAILED...
   fi
done


